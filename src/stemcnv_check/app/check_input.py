import importlib.resources
import os
import re
import ruamel.yaml as ruamel_yaml

from collections import defaultdict
from collections.abc import MutableMapping
from loguru import logger as logging

from stemcnv_check import STEM_CNV_CHECK, helpers
from stemcnv_check.exceptions import SampleFormattingError, SampleConstraintError, ConfigValueError


def flatten(dictionary, parent_key='', separator=':'):
    items = []
    for key, value in dictionary.items():
        new_key = parent_key + separator + key if parent_key else key
        if isinstance(value, MutableMapping):
            items.extend(flatten(value, new_key, separator=separator).items())
        else:
            # if isinstance(value, list):
            #     value = tuple(value)
            items.append((new_key, value))
    return dict(items)

# FIXME: scheme validation of hotspot tables!
# > maybe best to do via snakemake json config?


### Sanity checks ###
@logging.catch((FileNotFoundError, SampleConstraintError, SampleFormattingError, ConfigValueError), reraise=True)
def check_sample_table(args):

    sample_table_file = args.sample_table
    config_file = args.config

    if not os.path.isfile(sample_table_file):
        raise FileNotFoundError(f"Sample table file '{sample_table_file}' does not exist.")
    if not os.path.isfile(config_file):
        raise FileNotFoundError(f"Config file '{config_file}' does not exist.")

    sample_data_df = helpers.read_sample_table(sample_table_file, args.column_remove_regex)
    yaml = ruamel_yaml.YAML(typ='safe')
    with importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files').joinpath('default_config.yaml').open() as f:
        default_config = yaml.load(f)

    # Check sample_ids are unique
    non_unique_ids = sample_data_df['Sample_ID'][sample_data_df['Sample_ID'].duplicated()].to_list()
    if non_unique_ids:
        raise SampleConstraintError('The following Sample_IDs occur more than once: ' + ', '.join(non_unique_ids))

    # Check sex values exist & are formatted correctly
    samples = sample_data_df['Sex'].to_dict()
    if any(not s for s in samples.values()):
        missing = [sid for sid, s in samples.items() if not s]
        raise SampleConstraintError("Missing values for 'Sex' in the samplesheet. Affected samples: " + ', '.join(missing))
    elif not all(s in ('m', 'f') for s in map(lambda x: x[0].lower(), samples.values())):
        missing = [f"{sid}: {s}" for sid, s in samples.items() if not s[0].lower() in ('m', 'f')]
        raise SampleFormattingError("Not all values of the 'Sex' column in the samplesheet can be coerced to 'm' or 'f'. Affected samples: " + ', '.join(missing))

    # Check that all reference samples exist
    ref_samples = sample_data_df.set_index('Reference_Sample')['Sex'].to_dict()
    if '' in ref_samples:
        del ref_samples['']
    missing_refs = [ref for ref in ref_samples.keys() if ref not in samples.keys()]
    if missing_refs:
        raise SampleConstraintError("These 'Reference_Sample's do not also exist in the 'Sample_ID' column of the samplesheet: " + ', '.join(missing_refs))
    # Give warning if sex of reference and sample don't match
    sex_mismatch = [
        f"{s} ({sex})" for s, sex in samples.items() if
            sample_data_df.loc[s]['Reference_Sample'] and
            sex[0].lower() != samples[sample_data_df.loc[s]['Reference_Sample']][0].lower()
    ]
    if sex_mismatch:
        logging.error("The following samples have a different sex annotation than their Reference_Sample: " + ', '.join(sex_mismatch))
        raise SampleConstraintError("The following samples have a different sex annotation than their Reference_Sample: " + ', '.join(sex_mismatch))

    # This includes the global array definition config
    config = helpers.load_config(args, inbuilt_defaults=False)
    # check that the array_definition block is loaded
    if 'array_definition' not in config:
        raise ConfigValueError(
            "The 'array_definition' block is missing from the config file" +
            "." if args.no_cache else " and no arrays are defined in the globally."
        )

    # Check that Chip_Name & Chip_Pos match the sentrix wildcard regex (& aren't empty!)
    sample_data = helpers.read_sample_table(sample_table_file, args.column_remove_regex, return_type='list')
    for constraint, val in (('sample_id', 'sid'), ('sentrix_name', 'n'), ('sentrix_pos', 'p')):
        pattern = helpers.config_extract(['wildcard_constraints', constraint], config, default_config)
        mismatch = [sid for sid, n, p, _, _, _, _, _ in sample_data if not re.match('^' + pattern + '$', eval(val)) and eval(val)]
        if mismatch:
            raise SampleConstraintError(f"The '{constraint}' values for these samples do not fit the expected constraints: " + ', '.join(mismatch))

    # Check that all Array_Name values are defined in config, the example files will NOT pass this check
    array_names = set(sample_data_df['Array_Name'])
    array_names_not_config = [array for array in array_names if array not in config['array_definition']]
    if array_names_not_config:
        # logging.error("The following Array_Name values are not defined in the config file: " + ', '.join(array_names_not_config))
        raise SampleConstraintError("The following Array_Name values are not defined in the config file: " + ', '.join(array_names_not_config))

    # Check that all Chip_Name values have the same Array_Name
    n_arrays = sample_data_df.groupby(['Chip_Name']).nunique()['Array_Name']
    chip_name_nonuniq_array = list(n_arrays[n_arrays > 1].to_dict().keys())
    if chip_name_nonuniq_array:
        logging.error("The following Chip_Name values are listed with different Array_Names: " + ', '.join(chip_name_nonuniq_array))
        raise SampleConstraintError("The following Chip_Name values are listed with different Array_Names: " + ', '.join(chip_name_nonuniq_array))

    # # Check optional 'Regions of Interest' column, non-matching regions will be taken as gene_names
    # if 'Regions_of_Intertest' in sample_data_df.columns:
    #     for sample in samples:
    #         regions = sample_data_df.loc[sample]['Regions_of_Intertest'].split(';')
    #         checks = [re.match('^(.*|)?(chr)?[0-9XY]{1,2}:[0-9]+-[0-9]+$', region) for region in regions]
    #         if not all(checks):
    #             raise SampleFormattingError(f"The 'Region_of_Interest' entry for this sample is not properly formatted: {sample_data_df.loc[sample]['Regions_of_Intertest']['Sample_ID']}. Format: (NAME_)?(chr)?[CHR]:[start]-[end], separate multiple regions with only ';'.")


@logging.catch((FileNotFoundError, ConfigValueError), reraise=True)
def check_config(args, minimal_entries_only=False):

    config_file = args.config
    sample_table_file = args.sample_table

    if not os.path.isfile(config_file):
        raise FileNotFoundError(f"Config file '{config_file}' does not exist.")
    if not os.path.isfile(sample_table_file):
        raise FileNotFoundError(f"Sample table file '{sample_table_file}' does not exist.")

    config = helpers.load_config(args, inbuilt_defaults=False)
    config_with_defaults = helpers.load_config(args)
    cache_path = helpers.get_cache_dir(args, config_with_defaults)
    yaml = ruamel_yaml.YAML(typ='safe')
    with importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files').joinpath('default_config.yaml').open() as f:
        default_config = yaml.load(f)
    with importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files').joinpath('allowedvalues_config.yaml').open() as f:
        allowed_values = yaml.load(f)

    ## Check required entries

    # Files: static-data/*
    for array in config['array_definition'].keys():
        for req in allowed_values['array_definition']['__array'].keys():
            if req in ('penncnv_pfb_file', 'penncnv_GCmodel_file', 'array_density_file', 'array_gaps_file'):
                infostr = "\nYou can create it by running `stemcnv-check make-staticdata`"
                minimal_file = False
            else:
                infostr = ""
                minimal_file = True

            # csv is optional, skip if not defined or empty
            if req == 'csv_manifest_file' and not config['array_definition'][array].get(req):
                logging.warning(
                    f"Optional config entry is not used: array_definition:{array}:{req}. "
                    f"This will cause some probes to be unusable."
                )
                continue
            # All other array definition fields need to be defined
            if req not in config['array_definition'][array] or not config['array_definition'][array][req]:
                raise ConfigValueError(f"Required config entry is missing: array_definition:{array}:{req}")
            req_value = config['array_definition'][array][req]
            # Check basic file path issues and genome version 
            if req_value.startswith('~'):
                raise ConfigValueError(
                    "File paths in the config may not start with '~'. "
                    "Use absolute paths for anything outside the current directory."
                )
            if req == 'genome_version':
                if req_value not in ('hg38', 'hg19', 'GRCh38', 'GRCh37'):
                    raise ConfigValueError(
                        f"Genome version '{req_value}' for array '{array}' is not supported. Use 'hg38' or 'hg19'."
                    )
            # Check if file(s) exists, non-minimal/auto-generated files may optionally be missing
            elif (
                    minimal_file and not minimal_entries_only and
                    req_value == '__cache-default__'
            ):
                file = helpers.get_array_file(req, array, config, cache_path)
                if not os.path.isfile(file):
                    raise FileNotFoundError(
                        f"Array definition file from cache '{array}:{req}' does not exist: {file}."
                    )
            elif (
                    (minimal_file if minimal_entries_only else req_value != '__cache-default__') and
                    not os.path.isfile(req_value)
            ):
                raise FileNotFoundError(
                    f"Array definition file for '{array}:{req}' does not exist: {req_value}." + infostr
                )
        # Check that genome versions match Illumina syntax for manifest files
        # Can not be sure, this is a hard rule so only raise a warning
        expected_pattern = 'A1\\.(csv(.gz)?|bpm)$' \
            if config['array_definition'][array]['genome_version'] in ('hg19', 'GRCh37') \
            else 'A2\\.(csv(.gz)?|bpm)$'
        for manifest in ('csv_manifest_file', 'bpm_manifest_file'):
            if not re.search(expected_pattern, os.path.basename(config['array_definition'][array][manifest])):
                genome_version = config['array_definition'][array]['genome_version']
                logging.warning(
                    f"Manifest file '{config['array_definition'][array][manifest]}' does not match the expected pattern "
                    f"for '{genome_version}' genome.\n"
                    f"Illumina manifest files for {genome_version} should end in '{expected_pattern[0:2]}'."
                )

    # Folders: log, data, raw-input
    # and other settings
    for req in ('data_path', 'log_path', 'raw_data_folder'):
        try:
            if not config[req]:
                raise ConfigValueError(f"Required config entry is missing: {req}")
            if not os.path.isdir(config[req]):
                warn_str = f"Entry for required setting '{req}' is not an existing folder ({config[req]})! Creating it now."
                logging.warning(warn_str)
                os.makedirs(config[req], exist_ok=False)
        except KeyError:
            raise ConfigValueError(f"Required config entry is missing: {req}")

    if minimal_entries_only:
        return None

    # Other required settings: reports/*/filetype
    if 'reports' in config:
        for rep in config['reports']:
            if rep == '_default_':
                continue
            try:
                if not config['reports'][rep]['file_type']:
                    raise ConfigValueError(f"Required config entry is missing: reports:{rep}:file_type")
            except KeyError:
                raise ConfigValueError(f"Required config entry is missing: reports:{rep}:file_type")

    ## Full check of all other entries, specifically check value ranges
    def parse_scientific(n):
        if isinstance(n, (int, float)):
            return n
        if re.match('^[0-9.]+e[0-9]+$', n):
            n, e = n.split('e')
            out = float(n) * 10**int(e)
            return int(out) if int(out) == out else out
        else:
            ValueError(f"{n} can not be coerced to a number")

    with importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files').joinpath('label_name_definitions.yaml').open() as f:
        defined_label_names = yaml.load(f)
    defined_filtersets = set(helpers.config_extract(['settings', 'probe_filter_sets'], config, default_config).keys())
    defined_CNV_categories = set(
        helpers.config_extract(['evaluation_settings', 'CNV_call_labels'], config, default_config).keys()
    )
    defined_sampletable_cols = helpers.read_sample_table(sample_table_file, args.column_remove_regex).columns

    check_functions = defaultdict(lambda: lambda x, v: bool(re.match('^'+str(v)+'$', str(x))))
    check_functions['labels'] = lambda x, v: x in defined_label_names[v]

    def_functions = {
        'list': lambda x, v: type(x) == list,
        'str': lambda x, v: type(x) == str,
        'noneorstr': lambda x, v: x is None or type(x) == str,
        'int': lambda x, v: type(parse_scientific(x)) == int,
        'noneorfloat': lambda x, v: x is None or type(parse_scientific(x)) in (int, float),
        'float': lambda x, v: type(parse_scientific(x)) in (int, float),
        'bool': lambda x, v: type(x) == bool,
        'len': lambda x, v: len(x) == v,
        'le': lambda x, v: parse_scientific(x) <= v,
        'ge': lambda x, v: parse_scientific(x) >= v,
        'filterset': lambda x, v: x in defined_filtersets or x == '_default_',
        'filtersetnodefault': lambda x, v: x in defined_filtersets,
        'filtersetplusnone': lambda x, v: x in defined_filtersets or x == '_default_' or x == 'none',
        'cnvcallcategories': lambda x, v: x in defined_CNV_categories,
        'insamplesheet': lambda x, v: re.sub('^_+', '', x) in defined_sampletable_cols,
        'sectionsorall': lambda x, v: x == '__all__' or (len(x) > 1 and all(i in defined_label_names['report_sections'] for i in x)),
    }
    check_functions.update(def_functions)

    #Helper function:
    def formatfunc(s):
        strpart = s.rstrip('0123456789')
        if strpart in def_functions:
            numpart = int(s[len(strpart):]) if s[len(strpart):] else None
        elif strpart.startswith('labels:'):
            strpart = 'labels'
            numpart = s[len('labels:'):]
        else: # regex function
            numpart = strpart
        return strpart, numpart

    # Check all config entries
    errors = defaultdict(list)
    allowed_plotsections = defined_label_names['report_plotsections'] + ['_default_']
    for flatkey, config_value in flatten(config).items():
        # Need to change the key for variable config sections
        flatkey_ = re.sub('array_definition:[^:]+', 'array_definition:__array', flatkey)
        flatkey_ = re.sub('reports:[^:]+', 'reports:__report', flatkey_)
        flatkey_ = re.sub('tools:[^:]+', 'tools:__tool', flatkey_)
        flatkey_ = re.sub('evaluation_settings:CNV_call_labels:[^:]+', 'evaluation_settings:CNV_call_labels:__category', flatkey_)
        flatkey_ = re.sub('settings:probe_filter_sets:[^:]+', 'settings:probe_filter_sets:__filterset', flatkey_)
        flatkey_ = re.sub('reports:__report:call.data.and.plots:({})'.format('|'.join(allowed_plotsections)),
                          'reports:__report:call.data.and.plots:__plotsection', flatkey_)
        funcs = helpers.config_extract(flatkey_.split(':'), allowed_values, allowed_values)
        if funcs is None:
            # Happens for unknown config entries; helpers.config_extract gives a warning in this case
            #print(flatkey, config_value)
            continue
        funcs, *list_funcs = funcs.split('__')
        for func_key in funcs.split('_'):
            func_key, func_value = formatfunc(func_key)
            if not check_functions[func_key](config_value, func_value):
                #import pdb; pdb.set_trace()
                if isinstance(config_value, list):
                    config_value = tuple(config_value)
                errors[(flatkey, config_value)].append((func_key, func_value))
        for func_key in list_funcs:
            func_key, func_value = formatfunc(func_key)
            check = all(check_functions[func_key](i, func_value) for i in config_value)
            if not check:
                if isinstance(config_value, list):
                    config_value = tuple(config_value)
                errors[(flatkey, config_value)].append((func_key, func_value))

    help_strings = defaultdict(lambda: lambda v: "matching this regex: " + v)
    help_strings['labels'] = lambda v: f"one of the defined labels for '{v}': " + ', '.join(defined_label_names[v])
    help_strings.update({
        'list': lambda v: 'in a list',
        'str': lambda v: 'characters (add quotes for numbers or similar)',
        'int': lambda v: 'integers (whole numbers)',
        'noneorint': lambda v: 'integers (whole numbers) or null for an emtpy value',
        'float': lambda v: 'numbers',
        'bool': lambda v: 'booleans (True/False)',
        'len': lambda v: f"exactly {v} entries",
        'le': lambda v: f"values <={v}",
        'ge': lambda v: f"values >={v}",
        'filterset': lambda v: "the defined filtersets ({}) or _default_".format(', '.join(defined_filtersets)),
        'filtersetnodefault': lambda v: "only the defined filtersets ({}) but not '_default_'".format(', '.join(defined_filtersets)),
        'filtersetplusnone': lambda v: "the defined filtersets ({}), '_default_' or 'none'".format(', '.join(defined_filtersets)),
        'cnvcallcategories': lambda v: "the defined CNV call categories: " + ', '.join(defined_CNV_categories),
        'sectionsorall': lambda v: "either '__all__' or a list with of defined report sections: " + ', '.join(defined_label_names['allowed_sections']),
        'insamplesheet': lambda v: "column names of the samplesheet: " + ', '.join(defined_sampletable_cols)
    })
    if errors:
        for (flatkey, config_value), func_list in errors.items():
            warn_str = f"The config entry '{config_value}' for '{flatkey}' is invalid. Value(s) need to be "
            warn_str += ', and '.join([help_strings[func](val) for func, val in func_list]) + '.'
            logging.error(warn_str)
        raise ConfigValueError('The config contains values that are not allowed: ' + warn_str)

