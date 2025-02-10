import pytest
import importlib.resources
import textwrap
import ruamel.yaml as ruamel_yaml
from copy import deepcopy
from unittest.mock import MagicMock
from stemcnv_check.app.check_input import check_config, check_sample_table
from stemcnv_check import STEM_CNV_CHECK
from stemcnv_check.exceptions import ConfigValueError, SampleFormattingError, SampleConstraintError


@pytest.fixture
def minimal_config_block():
    return {
        'array_definition': {
            'default': {
                'genome_version': 'hg38',
                'bpm_manifest_file': 'static/manifest_A2.bpm',
                'egt_cluster_file': 'static/cluster.egt',
                'csv_manifest_file': '',
                'penncnv_pfb_file': 'dummy',
                'penncnv_GCmodel_file': 'dummy',
                'array_density_file': 'dummy',
                'array_gaps_file': 'dummy',
            }
        },
        'raw_data_folder': 'rawdata',
        'data_path': 'data',
        'log_path': 'logs',
    }

@pytest.fixture
def full_config_block(minimal_config_block):
    block = deepcopy(minimal_config_block)
    block['array_definition']['default'].update({
        'csv_manifest_file': 'static/manifest_A2.csv',
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
    if not fs.isdir('static'):
        fs.create_dir('static')
    for k, file in update_block['array_definition']['default'].items():
        if not file or file == 'genome_version' or file == 'dummy':
            continue
        if not fs.isfile(file):
            fs.create_file(file)


def test_check_sample_table(minimal_config_block, fs):

    default_config = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
    fs.add_real_file(default_config, read_only=True)

    base_sample_table = textwrap.dedent("""\
        Array_Name\tSample_ID\tChip_Name\tChip_Pos\tSex\tReference_Sample\tRegions_of_Interest\tSample_Group
        default\tSample1\t123456789000\tR01C01\tM\t\t\t\n
        default\tSample2\t123456789000\tR01C03\tM\tSample1\t\t\n
        """)

    testconfig = deepcopy(minimal_config_block)
    yaml = ruamel_yaml.YAML()
    with open('config.yaml', 'w') as f:
        yaml.dump(testconfig, f)
    sampletable = fs.create_file('sample_table.tsv', contents=base_sample_table)

    args = MagicMock(sample_table='sample_table.tsv', config='config.yaml', column_remove_regex=None)

    # Checks:
    # - success
    check_sample_table(args)

    # - dup ID
    sampletable.set_contents(base_sample_table + 'default\tSample1\t123456789000\tR01C01\tM\t\t\t\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table(args)

    # - missing sex
    sampletable.set_contents(base_sample_table + 'default\tSample3\t123456789000\tR01C03\t\tSample1\t\t\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table(args)

    # - sex mis-format
    sampletable.set_contents(base_sample_table + 'default\tSample3\t123456789000\tR01C03\twoman\tSample1\t\t\n')
    with pytest.raises(SampleFormattingError):
        check_sample_table(args)

    # - unknown ref
    sampletable.set_contents(base_sample_table + 'default\tSample3\t123456789000\tR01C03\tf\tSample0\t\t\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table(args)

    # - ref sex mismatch
    sampletable.set_contents(base_sample_table + 'default\tSample3\t123456789000\tR01C03\tf\tSample1\t\t\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table(args)

    # - wildcard constraint (mis)match
    sampletable.set_contents(base_sample_table + 'default\tSample3\tChipName\tR01C03\tm\tSample1\t\t\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table(args)

    # - array name not in config
    sampletable.set_contents(base_sample_table + 'undefined\tSample3\tChipName\tR01C03\tm\tSample1\t\t\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table(args)

    # - different array names for one chip_name
    testconfig['array_definition']['default2'] = dict()
    with open('config.yaml', 'w') as f:
        yaml.dump(testconfig, f)
    sampletable.set_contents(base_sample_table + 'default2\tSample3\t123456789000\tR01C05\tm\tSample1\t\t\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table(args)

    # Check for missing file
    args.sample_table = 'non_existent_sample_table.txt'
    with pytest.raises(FileNotFoundError):
        check_sample_table(args)


def test_check_config(minimal_config_block, full_config_block, fs, caplog):

    allowed_values = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'allowedvalues_config.yaml')
    default_config = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
    defined_labels = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'label_name_definitions.yaml')
    fs.add_real_file(default_config, read_only=True)
    fs.add_real_file(allowed_values, read_only=True)
    fs.add_real_file(defined_labels, read_only=True)
    fs.create_file(
        'sample_table.tsv',
        contents='Array_Name\tSample_ID\tChip_Name\tChip_Pos\tSex\tReference_Sample\tRegions_of_Interest\tSample_Group\n'
                 'default\nSample1\tChip1\t1\tM\tSample2\t\t\n'
    )
    fake_args = MagicMock(config='config.yaml', sample_table='sample_table.tsv', column_remove_regex=None)
    yaml = ruamel_yaml.YAML()
    testconfig = deepcopy(minimal_config_block)
    def update_config(testconfig):
        with open('config.yaml', 'w') as f:
            yaml.dump(testconfig, f)
    update_config(testconfig)

    # Error on missing config or sampletable file
    fake_args.config = 'non_existent_config.yaml'
    with pytest.raises(FileNotFoundError) as error:
        check_config(fake_args)
    assert 'non_existent_config.yaml' in str(error.value)

    fake_args.config = 'config.yaml'
    fake_args.sample_table = 'non_existent_sample_table.tsv'
    with pytest.raises(FileNotFoundError) as error:
        check_config(fake_args)
    assert 'non_existent_sample_table.tsv' in str(error.value)

    fake_args.sample_table = 'sample_table.tsv'

    # Check that the minimal config block works as long as only the required entries are checked
    # > bpm & egt files need to exist, other entries need to exist (except csv which can always be empty)
    prepare_fakefs_config(default_config, minimal_config_block, fs)
    assert not fs.isdir('data') and not fs.isdir('logs') and not fs.isdir('rawdata')
    check_config(fake_args, minimal_entries_only=True)
    # Also check that warnings are raised for missing folders and that they are created
    logrecords = caplog.records[-3:]
    assert [rec.levelname for rec in logrecords] == ['WARNING'] * 3
    assert [rec.message for rec in logrecords] == [
        f"Entry for required setting '{dirname}' is not an existing folder ({testconfig[dirname]})! Creating it now."
        for dirname in ('data_path', 'log_path', 'raw_data_folder')
    ]
    assert fs.isdir('data') and fs.isdir('logs') and fs.isdir('rawdata')

    # Check for fail on missing required entries
    del testconfig['array_definition']['default']['egt_cluster_file']
    update_config(testconfig)
    with pytest.raises(ConfigValueError):
        check_config(fake_args, minimal_entries_only=True)

    # Check for fail on missing required file
    testconfig['array_definition']['default']['egt_cluster_file'] = 'missing.egt'
    update_config(testconfig)
    with pytest.raises(FileNotFoundError) as error:
        check_config(fake_args)
    assert 'missing.egt' in str(error.value) and 'egt_cluster_file' in str(error.value)

    # csv file needs to exist even with minimal_entries_only, if it is defined
    testconfig = deepcopy(minimal_config_block)
    testconfig['array_definition']['default']['csv_manifest_file'] = 'static/manifest_A2.csv'
    update_config(testconfig)
    with pytest.raises(FileNotFoundError)as error:
        check_config(fake_args, minimal_entries_only=True)
    assert 'manifest_A2.csv' in str(error.value) and 'csv_manifest_file' in str(error.value)

    # check that warning for A1/A2 <> genome version is raised
    testconfig['array_definition']['default']['csv_manifest_file'] = 'static/manifest_A1.csv'
    fs.create_file('static/manifest_A1.csv')
    update_config(testconfig)
    check_config(fake_args, minimal_entries_only=True)
    logrecord = caplog.records[-1]
    assert logrecord.levelname == 'WARNING'
    assert logrecord.message == (
        "Manifest file 'static/manifest_A1.csv' does not match the expected pattern "
        "for 'hg38' genome.\nIllumina manifest files for hg38 should end in 'A2'."
    )

    # Check that without minimal_entries_only all static files need to exist
    testconfig = deepcopy(minimal_config_block)
    update_config(testconfig)
    with pytest.raises(FileNotFoundError) as error:
        check_config(fake_args)
    assert 'dummy' in str(error.value)

    # With all files existing, the check should pass
    prepare_fakefs_config(default_config, full_config_block, fs)
    check_config(fake_args)

    # Check for warning on unknown entries
    testconfig = deepcopy(full_config_block)
    testconfig['unknown_entry'] = 'nonsense'
    update_config(testconfig)
    check_config(fake_args)
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
        check_config(fake_args)
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
    # - labels:... functions
    # - Special other fields: filterset, filtersetdefault, sectionorsall
    # - check that all the smaller functions (len, ge, le, ...) work as expected


def test_default_config(full_config_block, caplog, fs):
    """Check that the actual default config passes the check_config function"""
    default_config = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
    allowed_values = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'allowedvalues_config.yaml')
    defined_labels = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'label_name_definitions.yaml')
    sample_table = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'sample_table_example.tsv')
    fs.add_real_file(default_config, read_only=True)
    fs.add_real_file(allowed_values, read_only=True)
    fs.add_real_file(defined_labels, read_only=True)
    fs.add_real_file(sample_table, read_only=True)

    prepare_fakefs_config(default_config, full_config_block, fs)
    fake_args = MagicMock(config='config.yaml', sample_table=sample_table, column_remove_regex=None)

    #Ensure that only expected warnings on missing folders are raised
    check_config(fake_args)
    logrecords = caplog.records
    assert len(logrecords) == 3
    assert [rec.message for rec in logrecords] == [
        "Entry for required setting 'data_path' is not an existing folder (data)! Creating it now.",
        "Entry for required setting 'log_path' is not an existing folder (logs)! Creating it now.",
        "Entry for required setting 'raw_data_folder' is not an existing folder (rawdata)! Creating it now."
    ]
