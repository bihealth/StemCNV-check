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
        'static_data': {
            'bpm_manifest_file': 'static/manifest.bpm',
            'egt_cluster_file': 'static/cluster.egt',
        },
        'raw_data_folder': 'rawdata',
        'data_path': 'data',
        'log_path': 'logs',
    }

@pytest.fixture
def full_config_block(minimal_config_block):
    block = minimal_config_block
    block['static_data'].update({
        'csv_manifest_file': 'static/manifest.csv',
        'genome_fasta_file': 'static/genome.fa',
        'genome_gtf_file': 'static/gencode.v45.gtf',
        'penncnv_pfb_file': 'static/chip.pfb',
        'penncnv_GCmodel_file': 'static/chip.gcmodel',
        'array_density_file': 'static/density.bed',
        'array_gaps_file': 'static/gaps.bed',
        'genomeInfo_file': 'static/genome_info.tsv',
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
    for k, file in update_block['static_data'].items():
        fs.create_file(file)


def test_check_sample_table(fs):

    default_config = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
    fs.add_real_file(default_config, read_only=True)

    base_sample_table = textwrap.dedent("""\
        Sample_ID\tChip_Name\tChip_Pos\tSex\tReference_Sample
        Sample1\t123456789000\tR01C01\tM\t\n
        Sample2\t123456789000\tR01C03\tM\tSample1\n
        """)

    fs.create_file('config.yaml', contents='data_path: data\nraw_data_folder: rawdata\nlog_path: logs\n')
    sampletable = fs.create_file('sample_table.tsv', contents=base_sample_table)

    # Check for missing file
    with pytest.raises(FileNotFoundError):
        check_sample_table('non_existent_sample_table.txt', 'non_existent_config.yaml')

    # Checks:
    # - success
    check_sample_table('sample_table.tsv', 'config.yaml')

    # - dup ID
    sampletable.set_contents(base_sample_table + 'Sample1\t123456789000\tR01C01\tM\t\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table('sample_table.tsv', 'config.yaml')

    # - missing sex
    sampletable.set_contents(base_sample_table + 'Sample3\t123456789000\tR01C03\t\tSample1\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table('sample_table.tsv', 'config.yaml')

    # - sex mis-format
    sampletable.set_contents(base_sample_table + 'Sample3\t123456789000\tR01C03\twoman\tSample1\n')
    with pytest.raises(SampleFormattingError):
        check_sample_table('sample_table.tsv', 'config.yaml')

    # - unknown ref
    sampletable.set_contents(base_sample_table + 'Sample3\t123456789000\tR01C03\tf\tSample0\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table('sample_table.tsv', 'config.yaml')

    # - ref sex mismatch
    sampletable.set_contents(base_sample_table + 'Sample3\t123456789000\tR01C03\tf\tSample1\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table('sample_table.tsv', 'config.yaml')

    # - wildcard constraint (mis)match
    sampletable.set_contents(base_sample_table + 'Sample3\tChipName\tR01C03\tm\tSample1\n')
    with pytest.raises(SampleConstraintError):
        check_sample_table('sample_table.tsv', 'config.yaml')

    #FIXME: - ROI formatting (needs to be adapted to new format allowing gband, gene or pos-string)



def test_check_config(minimal_config_block, full_config_block, fs, caplog):

    allowed_values = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'allowedvalues_config.yaml')
    default_config = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
    fs.add_real_file(default_config, read_only=True)
    fs.add_real_file(allowed_values, read_only=True)
    fs.create_file('sample_table.tsv',
                   contents='Sample_ID\tChip_Name\tChip_Pos\tSex\tReference_Sample\nSample1\tChip1\t1\tM\tSample2\n')

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
    testconfig['static_data']['egt_cluster_file'] = 'missing.bpm'
    update_config(testconfig)
    with pytest.raises(FileNotFoundError):
        check_config('config.yaml', 'sample_table.tsv')

    # Check for fail on missing required entries
    del testconfig['static_data']['egt_cluster_file']
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
    # - wrong type
    testconfig['settings'] = {'chromosomes': '1-22'}
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
        "The config entry '1-22' for 'settings:chromosomes' is invalid. " +
        "Value(s) need to be in a list, and matching this regex: (chr)?[0-9XY]+.",
        "The config entry '1.5' for 'settings:array_attribute_summary:density.windows' is invalid. " +
        "Value(s) need to be integers (whole numbers).",
        "The config entry '('PennCNV', 'GATK')' for 'settings:CNV.calling.tools' is invalid. " +
        "Value(s) need to be matching this regex: (PennCNV|CBS).",
    ]

    #FIXME: future tests
    # - Check for matching/mismtahcing sample_table columns
    # - Special other fields: filterset, filtersetdefault, section, sectionsall
    # - check that all the smaller functions (len, ge, le, ...) work as expected


def test_default_config(full_config_block, fs):
    """Check that the actual default cnfig passes the check_config function"""
    default_config = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
    allowed_values = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'allowedvalues_config.yaml')
    sample_table = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'sample_table_example.tsv')
    fs.add_real_file(default_config, read_only=True)
    fs.add_real_file(allowed_values, read_only=True)
    fs.add_real_file(sample_table, read_only=True)

    prepare_fakefs_config(default_config, full_config_block, fs)

    fs.listdir('static')
    check_config('config.yaml', sample_table)
