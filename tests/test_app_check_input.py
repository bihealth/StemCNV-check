import pytest
import importlib.resources
import textwrap
import ruamel.yaml as ruamel_yaml
from copy import deepcopy
from stemcnv_check.app.check_input import check_config, check_sample_table
from stemcnv_check import STEM_CNV_CHECK
from stemcnv_check.exceptions import ConfigValueError, SampleFormattingError, SampleConstraintError


@pytest.fixture
def minimal_config_block():
    return {
        'array_definition': {
            'default': {
                'genome_version': 'hg38',
                'bpm_manifest_file': 'static/manifest.bpm',
                'egt_cluster_file': 'static/cluster.egt',
            }
        },
        'raw_data_folder': 'rawdata',
        'data_path': 'data',
        'log_path': 'logs',
    }

@pytest.fixture
def full_config_block(minimal_config_block):
    block = minimal_config_block
    block['array_definition']['default'].update({
        'csv_manifest_file': 'static/manifest.csv',
        'penncnv_pfb_file': 'static/chip.pfb',
        'penncnv_GCmodel_file': 'static/chip.gcmodel',
        'array_density_file': 'static/density.bed',
        'array_gaps_file': 'static/gaps.bed',
    })
    return block


def prepare_fakefs_config(default_config, update_block, fs):
    # Update required entries and write to fs / 'config.yaml'
    yaml = ruamel_yaml.YAML()
    with open(default_config, 'r') as f:
        default_config = yaml.load(f)
    default_config.update(update_block)
    with open('config.yaml', 'w') as f:
        yaml.dump(default_config, f)
    # Create the fake files:
    fs.create_dir('static')
    for k, file in update_block['array_definition']['default'].items():
        if file == 'genome_version':
            continue
        fs.create_file(file)


def test_check_sample_table(minimal_config_block, fs):

    default_config = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
    fs.add_real_file(default_config, read_only=True)

    base_sample_table = textwrap.dedent("""\
        Array_Name\tSample_ID\tChip_Name\tChip_Pos\tSex\tReference_Sample
        default\tSample1\t123456789000\tR01C01\tM\t\n
        default\tSample2\t123456789000\tR01C03\tM\tSample1\n
        """)

    testconfig = deepcopy(minimal_config_block)
    yaml = ruamel_yaml.YAML()
    with open('config.yaml', 'w') as f:
        yaml.dump(testconfig, f)
    sampletable = fs.create_file('sample_table.tsv', contents=base_sample_table)

    # Check for missing file
    with pytest.raises(FileNotFoundError):
        check_sample_table('non_existent_sample_table.txt', 'non_existent_config.yaml')

    # Checks:
    # - success
    check_sample_table('sample_table.tsv', 'config.yaml')

    # - dup ID
    sampletable.set_contents(base_sample_table + 'default\tSample1\t123456789000\tR01C01\tM\t\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table('sample_table.tsv', 'config.yaml')

    # - missing sex
    sampletable.set_contents(base_sample_table + 'default\tSample3\t123456789000\tR01C03\t\tSample1\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table('sample_table.tsv', 'config.yaml')

    # - sex mis-format
    sampletable.set_contents(base_sample_table + 'default\tSample3\t123456789000\tR01C03\twoman\tSample1\n')
    with pytest.raises(SampleFormattingError):
        check_sample_table('sample_table.tsv', 'config.yaml')

    # - unknown ref
    sampletable.set_contents(base_sample_table + 'default\tSample3\t123456789000\tR01C03\tf\tSample0\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table('sample_table.tsv', 'config.yaml')

    # - ref sex mismatch
    sampletable.set_contents(base_sample_table + 'default\tSample3\t123456789000\tR01C03\tf\tSample1\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table('sample_table.tsv', 'config.yaml')

    # - wildcard constraint (mis)match
    sampletable.set_contents(base_sample_table + 'default\tSample3\tChipName\tR01C03\tm\tSample1\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table('sample_table.tsv', 'config.yaml')

    # - array name not in config
    sampletable.set_contents(base_sample_table + 'undefined\tSample3\tChipName\tR01C03\tm\tSample1\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table('sample_table.tsv', 'config.yaml')

    # - different array names for one chip_name
    testconfig['array_definition']['default2'] = dict()
    with open('config.yaml', 'w') as f:
        yaml.dump(testconfig, f)
    sampletable.set_contents(base_sample_table + 'default2\tSample3\t123456789000\tR01C05\tm\tSample1\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table('sample_table.tsv', 'config.yaml')


def test_check_config(minimal_config_block, full_config_block, fs, caplog):

    allowed_values = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'allowedvalues_config.yaml')
    default_config = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
    fs.add_real_file(default_config, read_only=True)
    fs.add_real_file(allowed_values, read_only=True)
    fs.create_file('sample_table.tsv',
                   contents='Array_Name\tSample_ID\tChip_Name\tChip_Pos\tSex\tReference_Sample\ndefault\nSample1\tChip1\t1\tM\tSample2\n')

    yaml = ruamel_yaml.YAML()
    testconfig = deepcopy(minimal_config_block)
    def update_config(testconfig):
        with open('config.yaml', 'w') as f:
            yaml.dump(testconfig, f)

    # Error on missing config file
    with pytest.raises(FileNotFoundError):
        check_config('non_existent_config.yaml', 'non_existent_sample_table.txt')

    # Check minimal config block with req_only works (entries need to exist, but not the files)
    # Also check that warnings are raised for missing folders and that they are created
    assert not fs.isdir('data') and not fs.isdir('logs') and not fs.isdir('rawdata')
    update_config(testconfig)
    check_config('config.yaml', 'sample_table.tsv', required_only=True)
    logrecords = caplog.records[-3:]
    assert [rec.levelname for rec in logrecords] == ['WARNING'] * 3
    assert [rec.message for rec in logrecords] == [
        f"Entry for required setting '{dirname}' is not an existing folder ({testconfig[dirname]})! Creating it now."
        for dirname in ('data_path', 'log_path', 'raw_data_folder')
    ]
    assert fs.isdir('data') and fs.isdir('logs') and fs.isdir('rawdata')

    # Check for fail on missing required files
    testconfig['array_definition']['default']['egt_cluster_file'] = 'missing.bpm'
    update_config(testconfig)
    with pytest.raises(FileNotFoundError):
        check_config('config.yaml', 'sample_table.tsv')

    # Check for fail on missing required entries
    del testconfig['array_definition']['default']['egt_cluster_file']
    update_config(testconfig)
    with pytest.raises(ConfigValueError):
        check_config('config.yaml', 'sample_table.tsv', required_only=True)

    # Check without req_only all static files need to exist
    with pytest.raises(FileNotFoundError):
        check_config('config.yaml', 'sample_table.tsv')

    # With files existing, the check should pass
    prepare_fakefs_config(default_config, full_config_block, fs)
    check_config('config.yaml', 'sample_table.tsv')

    # Check for warning on unknown entries
    testconfig = deepcopy(full_config_block)
    testconfig['unknown_entry'] = 'nonsense'
    update_config(testconfig)
    check_config('config.yaml', 'sample_table.tsv')
    logrecord = caplog.records[-1]
    assert logrecord.levelname == 'WARNING'
    assert logrecord.message == '"unknown_entry" is not a valid config entry or has been deprecated'

    # Check for Error on entry outside specifications:
    del testconfig['unknown_entry']
    testconfig['settings'] = dict()
    # - wrong type
    testconfig['settings']['PennCNV'] = {'enable_LOH_calls': 'True'}
    # - wrong number type (float vs int)
    testconfig['settings']['array_attribute_summary'] = {'density.windows': 1.5}
    # - value not matching regex
    testconfig['settings']['CNV.calling.tools'] = ['PennCNV', 'GATK']
    update_config(testconfig)
    with pytest.raises(ConfigValueError):
        check_config('config.yaml', 'sample_table.tsv')
    logrecords = caplog.records[-3:]
    assert [rec.levelname for rec in logrecords] == ['ERROR'] * 3
    assert [rec.message for rec in logrecords] == [
        "The config entry 'True' for 'settings:PennCNV:enable_LOH_calls' is invalid. " +
        "Value(s) need to be booleans (True/False).",
        "The config entry '1.5' for 'settings:array_attribute_summary:density.windows' is invalid. " +
        "Value(s) need to be integers (whole numbers).",
        "The config entry '('PennCNV', 'GATK')' for 'settings:CNV.calling.tools' is invalid. " +
        "Value(s) need to be matching this regex: (PennCNV|CBS).",
    ]

    #FIXME: future tests
    # - Check for matching/mismatching sample_table columns
    # - Special other fields: filterset, filtersetdefault, section, sectionsall
    # - check that all the smaller functions (len, ge, le, ...) work as expected


def test_default_config(full_config_block, caplog, fs):
    """Check that the actual default config passes the check_config function"""
    default_config = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
    allowed_values = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'allowedvalues_config.yaml')
    sample_table = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'sample_table_example.tsv')
    fs.add_real_file(default_config, read_only=True)
    fs.add_real_file(allowed_values, read_only=True)
    fs.add_real_file(sample_table, read_only=True)

    prepare_fakefs_config(default_config, full_config_block, fs)

    #Ensure that only expected warnings on missing folders are raised
    check_config('config.yaml', sample_table)
    logrecords = caplog.records
    assert len(logrecords) == 3
    assert [rec.message for rec in logrecords] == [
        "Entry for required setting 'data_path' is not an existing folder (data)! Creating it now.",
        "Entry for required setting 'log_path' is not an existing folder (logs)! Creating it now.",
        "Entry for required setting 'raw_data_folder' is not an existing folder (rawdata)! Creating it now."
    ]
