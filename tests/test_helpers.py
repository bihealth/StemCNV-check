import pytest
import textwrap
import stemcnv_check.helpers as helpers
import importlib.resources
from stemcnv_check import STEM_CNV_CHECK
from stemcnv_check.exceptions import ConfigValueError, SampleConstraintError
from unittest.mock import patch, MagicMock
from collections import OrderedDict

@pytest.fixture
def sample_table_minimal():
    return textwrap.dedent(
        """\
        Chip_Name\tChip_Pos\tSample_ID\tSex\tReference_Sample
        123456789000\tR01C01\tCellline-A-MB\tFemale\t
        123456789001\tR01C01\tCellline-A-WB\tFemale\tCellline-A-MB
        """
    )

@pytest.fixture
def sample_table_missing():
    return textwrap.dedent(
        """\
        Chip_Name\tSample_ID\tSex\tReference_Sample
        123456789000\tCellline-A-MB\tFemale\t
        123456789001\tCellline-A-WB\tFemale\tCellline-A-MB
        """
    )

@pytest.fixture
def sample_table_extra_cols():
    return textwrap.dedent(
        """\
        # Commented line
        Sample_Name\tTest_col\tSex\tReference_Sample\tSample_Group\tRegions_of_Interest\tChip_Name\tChip_Pos\tSample_ID
        Cellline-A MasterBank\t123,456\tFemale\t\tExampleCellines\tExample1|chr1:100000-200000\t123456789000\tR01C01\tCellline-A-MB
        Cellline-A WorkingBank\tCellline-B-MB\tFemale\tCellline-A-MB\tExampleCellines\tExample1|chr1:100000-200000;chr11:60000-70000\t123456789001\tR01C01\tCellline-A-WB
        Cellline-B MasterBank\t0\tMale\t\tExampleCellines\tExample1|chr1:100000-200000;chr11:60000-70000\t123456789000\tR01C02\tCellline-B-MB
        #removed sample\t0\tMale\t\tExampleCellines\t\t123456789000\t000\tabc
        Cellline-B-1 clone1\t\tMale\tCellline-B-MB\tExampleCellines\t\t123456789001\tR01C02\tCellline-B-1-cl1
        """)


def test_read_sample_table(sample_table_minimal, sample_table_missing, sample_table_extra_cols, fs):
    expected = [["Cellline-A-MB", "123456789000", "R01C01", "Female", ""],
                ["Cellline-A-WB", "123456789001", "R01C01", "Female", "Cellline-A-MB"]]

    fs.create_file('sample_table_minimal.tsv', contents=sample_table_minimal)
    fs.create_file('sample_table_missing.tsv', contents=sample_table_missing)
    fs.create_file('sample_table_extra_cols.tsv', contents=sample_table_extra_cols)

    # test minimal sample table
    assert expected == helpers.read_sample_table("sample_table_minimal.tsv")

    # test minimal with dict output
    cols = ['Sample_ID', 'Chip_Name', 'Chip_Pos', 'Sex', 'Reference_Sample']
    expected_dict = [OrderedDict((k, v) for k, v in zip(cols, data_row)) for data_row in expected]
    assert expected_dict == helpers.read_sample_table("sample_table_minimal.tsv", with_opt=True)

    # test failer with missing columns
    with pytest.raises(SampleConstraintError):
        helpers.read_sample_table("sample_table_missing.tsv")

    # test extended example without dict output
    expected += [["Cellline-B-MB", "123456789000", "R01C02", "Male", ""],
                 ["Cellline-B-1-cl1", "123456789001", "R01C02", "Male", "Cellline-B-MB"]]
    assert expected == helpers.read_sample_table("sample_table_extra_cols.tsv")

    # test extended example without dict output
    extra_expected = [
        ["Cellline-A MasterBank", "123,456", "ExampleCellines", "Example1|chr1:100000-200000"],
        ["Cellline-A WorkingBank", "Cellline-B-MB", "ExampleCellines", "Example1|chr1:100000-200000;chr11:60000-70000"],
        ["Cellline-B MasterBank", "0", "ExampleCellines", "Example1|chr1:100000-200000;chr11:60000-70000"],
        ["Cellline-B-1 clone1", "", "ExampleCellines", ""]
    ]
    cols_extra = cols + ["Sample_Name", "Test_col", "Sample_Group", "Regions_of_Interest"]
    expected_dict = [OrderedDict((k, v) for k, v in zip(cols_extra, data_row+extra_row))
                     for data_row, extra_row in zip(expected, extra_expected)]
    assert expected_dict == helpers.read_sample_table("sample_table_extra_cols.tsv", with_opt=True)

    # Test the included example file
    fs.add_real_file(
        importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'sample_table_example.tsv'),
        read_only=True
    )
    cols_extra = cols + ["Sample_Name", "Sample_Group", "Regions_of_Interest"]
    expected_dict = [OrderedDict((k, v) for k, v in zip(cols_extra, data_row+[extra_row[0]]+extra_row[2:]))
                     for data_row, extra_row in zip(expected, extra_expected)]
    assert expected_dict == helpers.read_sample_table(
        importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'sample_table_example.tsv'),
        with_opt=True
    )



@patch('stemcnv_check.helpers.importlib.resources.files')
def test_make_singularity_args(mock_resource_files, fs):
    mock_resource_files.return_value = '/fake/snakedir'

    config = {'data_path': '/path/to/data',
              'raw_data_folder': '/path/to/rawdata',
              'log_path': '/path/to/logs',
              'static_data': {
                  'genome_fasta_file': 'relative/genome.fasta'
              }}
    expected = "-B /path/to/data:/outside/data,/path/to/rawdata:/outside/rawdata,/path/to/logs:/outside/logs," + \
               "/fake/snakedir:/outside/snakedir"
    assert expected == helpers.make_apptainer_args(config, not_existing_ok=True)

    with pytest.raises(FileNotFoundError):
        helpers.make_apptainer_args(config)

    assert expected + ',/tmp/tmpdir:/outside/tmp' == helpers.make_apptainer_args(config,
                                                                                 tmpdir='/tmp/tmpdir', not_existing_ok=True)


    fs.create_file('relative/genome.fasta')
    expected += ",relative/genome.fasta:/outside/static/genome.fasta"
    assert expected == helpers.make_apptainer_args(config)

@pytest.fixture
def user_config():
    return textwrap.dedent(
        """\
        data_path: /path/to/data
        cat:
            key: value
            nested:
                key: value2
        """
    )

@pytest.fixture
def default_config():
    return textwrap.dedent(
        """\
        data_path: /path/to/data
        cat:
            key: othervalue
            nested:
                key: value2
                key2: val3
        """
    )


def test_load_config(user_config, default_config, fs):

    fs.create_file('config.yaml', contents=user_config)
    fs.create_file(
        importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml'),
        contents=default_config
    )

    config_dict = {
        'data_path': '/path/to/data',
        'cat': {'key': 'value',
                'nested': {'key': 'value2'}}}
    config_from_defaults = {
        'data_path': '/path/to/data',
        'cat': {'key': 'value',
                'nested': {'key': 'value2',
                           'key2': 'val3'}}}

    assert config_dict == helpers.load_config('config.yaml', False)
    assert config_from_defaults == helpers.load_config('config.yaml', True)


def test_config_extract(caplog):
    config_dict   = {'data_path': '/path/to/data',
                     'cat': {'key': 'value',
                             'nested': {'key': 'value2'}}}
    default_config = {'data_path': '/path/to/data',
                      'cat': {'key': 'othervalue',
                              'nested': {'key': 'value2',
                                         'key2': 'val3'}}}

    assert '/path/to/data' == helpers.config_extract(('data_path',), config_dict, default_config)

    assert 'val3' == helpers.config_extract(('cat', 'nested', 'key2'), config_dict, default_config)
    assert caplog.records[-1].levelname == 'DEBUG'
    assert "Using config default values for: cat : nested : key2" == caplog.records[-1].message

    assert 'value' == helpers.config_extract(('cat', 'key'), config_dict, default_config)

    assert helpers.config_extract(('cat', 'nested', 'key3'), config_dict, default_config) is None
    assert caplog.records[-1].levelname == 'WARNING'
    assert '"cat : nested : key3" is not a valid config entry or has been deprecated' == caplog.records[-1].message


def test_collect_SNP_cluster_ids(sample_table_extra_cols, fs):
    fs.create_file('sample_table.tsv', contents=sample_table_extra_cols)
    sample_data_full = helpers.read_sample_table('sample_table.tsv', with_opt=True)

    all_ids = ['Cellline-A-MB', 'Cellline-A-WB', 'Cellline-B-MB', 'Cellline-B-1-cl1']

    # Test finding sample_ids based on machting entries in given column
    # The search sample ID itself is *never* included in the result
    assert set(all_ids[1:]) == helpers.collect_SNP_cluster_ids('Cellline-A-MB', ['__Sample_Group'], sample_data_full)
    assert {all_ids[0]} == helpers.collect_SNP_cluster_ids('Cellline-A-WB', ['__Chip_Pos'], sample_data_full)
    with pytest.raises(ConfigValueError):
        helpers.collect_SNP_cluster_ids('Cellline-A-MB', ['__NonExisting'], sample_data_full)
    # Test finding sample_ids based on values in given column
    assert {'123', '456'} == helpers.collect_SNP_cluster_ids('Cellline-A-MB', ['_Test_col'], sample_data_full)
    assert {all_ids[2]} == helpers.collect_SNP_cluster_ids('Cellline-A-WB', ['_Test_col'], sample_data_full)
    with pytest.raises(ConfigValueError):
        helpers.collect_SNP_cluster_ids('Cellline-A-WB', ['_Testcol'], sample_data_full)
    # Test using sample_ids directly
    assert {all_ids[2], 'Test'} == helpers.collect_SNP_cluster_ids('Cellline-A-WB', ['_Test_col', 'Test'], sample_data_full)
