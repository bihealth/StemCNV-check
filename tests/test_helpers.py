import os
import pytest
import textwrap
import stemcnv_check.helpers as helpers
import importlib.resources
from pathlib import Path
from stemcnv_check import STEM_CNV_CHECK
from stemcnv_check.exceptions import ConfigValueError, SampleConstraintError, CacheUnavailableError
from unittest.mock import patch, MagicMock
from collections import defaultdict, OrderedDict

@pytest.fixture
def sample_table_minimal():
    return textwrap.dedent(
        """\
        Array_Name\tChip_Name\tChip_Pos\tSample_ID\tSex\tReference_Sample\tRegions_of_Interest\tSample_Group
        ExampleArray\t123456789000\tR01C01\tCellline-A-MB\tFemale\t\tExample1|chr1:100000-200000\tExampleCellines
        ExampleArray\t123456789001\tR01C01\tCellline-A-WB\tFemale\tCellline-A-MB\tExample1|chr1:100000-200000;chr11:60000-70000\tExampleCellines
        """
    )

@pytest.fixture
def sample_table_missing():
    return textwrap.dedent(
        """\
        Array_Name\tChip_Name\tSample_ID\tSex\tReference_Sample\tRegions_of_Interest\tSample_Group
        ExampleArray\t123456789000\tCellline-A-MB\tFemale\t\tExample1|chr1:100000-200000\tExampleCellines
        ExampleArray\t123456789001\tCellline-A-WB\tFemale\tCellline-A-MB\tExample1|chr1:100000-200000;chr11:60000-70000\tExampleCellines
        """
    )

@pytest.fixture
def sample_table_extra_cols():
    return textwrap.dedent(
        """\
        # Commented line
        Array_Name\tSample_Name\tTest_col\tSex\tReference_Sample\tRegions_of_Interest\tSample_Group\tChip_Name\tChip_Pos\tSample_ID
        ExampleArray\tCellline-A MasterBank\tCellline-A-WB,Cellline-B-MB\tFemale\t\tExample1|chr1:100000-200000\tExampleCellines\t123456789000\tR01C01\tCellline-A-MB
        ExampleArray\tCellline-A WorkingBank\tCellline-B-MB\tFemale\tCellline-A-MB\tExample1|chr1:100000-200000;chr11:60000-70000\tExampleCellines\t123456789001\tR01C01\tCellline-A-WB
        ExampleArray\tCellline-B MasterBank\t0\tMale\t\tExample1|chr1:100000-200000;chr11:60000-70000\tExampleCellines\t123456789000\tR01C02\tCellline-B-MB
        #ExampleArray\tremoved sample\t0\tMale\t\t\tExampleCellines\t123456789000\t000\tabc
        ExampleArray\tCellline-B-1 clone1\t\tMale\tCellline-B-MB\t\tExampleCellines\t123456789001\tR01C02\tCellline-B-1-cl1
        """)


def test_read_sample_table(sample_table_minimal, sample_table_missing, sample_table_extra_cols, fs):
    expected = [
        ["Cellline-A-MB", "123456789000", "R01C01", "ExampleArray", "Female", "", "Example1|chr1:100000-200000", "ExampleCellines"],
        ["Cellline-A-WB", "123456789001", "R01C01", "ExampleArray", "Female", "Cellline-A-MB", "Example1|chr1:100000-200000;chr11:60000-70000", "ExampleCellines"]
    ]

    fs.create_file('sample_table_minimal.tsv', contents=sample_table_minimal)
    fs.create_file('sample_table_missing.tsv', contents=sample_table_missing)
    fs.create_file('sample_table_extra_cols.tsv', contents=sample_table_extra_cols)

    # test minimal sample table
    assert expected == helpers.read_sample_table("sample_table_minimal.tsv", return_type='list')

    # test minimal with dict output
    cols = [
        'Sample_ID', 'Chip_Name', 'Chip_Pos', 'Array_Name', 'Sex',
        'Reference_Sample', 'Regions_of_Interest', 'Sample_Group'
    ]
    expected_dict = [OrderedDict((k, v) for k, v in zip(cols, data_row)) for data_row in expected]
    assert expected_dict == helpers.read_sample_table("sample_table_minimal.tsv", return_type='list_withopt')

    # test fail with missing columns
    with pytest.raises(SampleConstraintError):
        helpers.read_sample_table("sample_table_missing.tsv")

    # test that column name removal with regex works
    assert expected == helpers.read_sample_table("sample_table_minimal.tsv", name_remove_regex=' .*', return_type='list')
    file_content = ('\t'.join([col + ' dummytext' for col in cols]) + '\n' +
                    '\n'.join(['\t'.join(line) for line in expected]) + '\n')
    fs.create_file('sample_table_colwithspaces.tsv', contents=file_content)
    assert expected == helpers.read_sample_table("sample_table_colwithspaces.tsv", name_remove_regex=' .*', return_type='list')

    # Test dataframe output
    import pandas as pd
    expected_df = pd.DataFrame(
        data=expected,
        columns=cols
    ).set_index('Sample_ID', drop=False)
    actual_df = helpers.read_sample_table("sample_table_minimal.tsv")
    assert expected_df.equals(actual_df)

    # test extended default without dict output
    expected += [
        ["Cellline-B-MB", "123456789000", "R01C02", "ExampleArray", "Male", "", "Example1|chr1:100000-200000;chr11:60000-70000", "ExampleCellines"],
        ["Cellline-B-1-cl1", "123456789001", "R01C02", "ExampleArray", "Male", "Cellline-B-MB", "", "ExampleCellines"]
    ]
    assert expected == helpers.read_sample_table("sample_table_extra_cols.tsv", return_type='list')

    # test extended example without dict output
    extra_expected = [
        ["Cellline-A MasterBank", "Cellline-A-WB,Cellline-B-MB"],
        ["Cellline-A WorkingBank", "Cellline-B-MB"],
        ["Cellline-B MasterBank", "0"],
        ["Cellline-B-1 clone1", ""]
    ]
    cols_extra = cols + ["Sample_Name", "Test_col"]
    expected_dict = [OrderedDict((k, v) for k, v in zip(cols_extra, data_row+extra_row))
                     for data_row, extra_row in zip(expected, extra_expected)]
    assert expected_dict == helpers.read_sample_table("sample_table_extra_cols.tsv", return_type='list_withopt')

    # Test that xlsx files work
    fs.add_real_file('tests/data/sample_table.xlsx', read_only=True)
    assert expected_dict == helpers.read_sample_table("tests/data/sample_table.xlsx", return_type='list_withopt')

    # Test the included example file
    fs.add_real_file(
        importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'sample_table_example.tsv'),
        read_only=True
    )
    expected_dict = [OrderedDict((k, v) for k, v in zip(cols, data_row+[extra_row[0]]+extra_row[2:]))
                     for data_row, extra_row in zip(expected, extra_expected)]
    assert expected_dict == helpers.read_sample_table(
        importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'sample_table_example.tsv'),
        return_type='list_withopt'
    )


@patch('stemcnv_check.helpers.importlib.resources.files')
def test_make_singularity_args(mock_resource_files, fs):
    mock_resource_files.return_value = '/fake/snakedir'
    cache_path = '/path/to/cache'
    config = {
        'data_path': '/path/to/data',
        'raw_data_folder': '/path/to/rawdata',
        'log_path': '/path/to/logs',
        'global_settings': {
            'hg19_genome_fasta': 'relative/genome.fasta',
            'hg19_gtf_file': 'relative/gtf.gtf',
            'hg19_genomeInfo_file': 'relative/genome.info',
            'hg19_mehari_transcript_db': 'relative/mehari_db',
            'dosage_sensitivity_scores': 'relative/dosage_sensitivity_scores.tsv',
            'cache_dir': cache_path
        },
        'array_definition': {
            'ExampleArray': {
                'genome_version': 'hg19',
                'bpm_manifest': 'relative_array/bpm_manifest.bpm',
                'penncnv_pfb_file': 'relative_array/pfb_file.pfb',
            }
        }
    }
    # 'make_singularity_args' expects to-be-bound paths to exist
    # output will be sorted by second path and include only absolute paths (with resolved links)
    # the fake fs puts everything in an empty root directory
    expected_base = [
        "'/path/to/data':'/outside/data'",
        "'/path/to/logs':'/outside/logs'",
        "'/path/to/rawdata':'/outside/rawdata'",
        "'/fake/snakedir':'/outside/snakedir'"
    ]
    def get_expected(extra):
        return "-B " + ','.join(sorted(expected_base + extra, key=lambda x: x.split(':')[1]))
    expected_extra = [
        "'/relative/genome.fasta':'/outside/static/genome.fasta'",
        "'/relative/genome.info':'/outside/static/genome.info'",
        "'/relative/gtf.gtf':'/outside/static/gtf.gtf'",
        "'/relative/mehari_db':'/outside/static/mehari_db'",
        "'/relative/dosage_sensitivity_scores.tsv':'/outside/static/dosage_sensitivity_scores.tsv'",
        "'/relative_array/bpm_manifest.bpm':'/outside/ExampleArray/bpm_manifest.bpm'"
    ]
    for file in expected_base + expected_extra[:-1]:
        fs.create_file(file.split(':')[0].rstrip("'").lstrip("'"))
    fs.create_dir('relative_array')

    # Fail if defined required file does not exist (bpm), even with not_existing_ok
    with pytest.raises(FileNotFoundError):
        helpers.make_apptainer_args(config, None, not_existing_ok=True)
    fs.create_file('/relative_array/bpm_manifest.bpm')

    # Successful test with not_existing_ok (pfb file missing)
    assert (get_expected(expected_extra+["'/relative_array':'/outside/writearray/ExampleArray'"]) ==
            helpers.make_apptainer_args(config, None, not_existing_ok=True))

    # Successful test with not_existing_ok & tmpdir (pfb file missing)
    fs.create_dir('/tmp/tmpdir')
    assert (get_expected(expected_extra + ["'/relative_array':'/outside/writearray/ExampleArray'", "'/tmp/tmpdir':'/outside/tmp'"])
            == helpers.make_apptainer_args(config, None, tmpdir='/tmp/tmpdir', not_existing_ok=True))

    # Fail if the defined file does not exist without not_existing_ok (pfb file)
    with pytest.raises(FileNotFoundError):
        helpers.make_apptainer_args(config, None)
    # test after file creation
    fs.create_file('relative_array/pfb_file.pfb')
    expected = get_expected(expected_extra + ["'/relative_array/pfb_file.pfb':'/outside/ExampleArray/pfb_file.pfb'"])
    assert expected == helpers.make_apptainer_args(config, None)

    # test with direct addition of mount path, these are NOT checked for existence
    assert expected+',/abcdef' == helpers.make_apptainer_args(config, None, extra_bind_args='/abcdef')

    # test that links are resolved to actual file paths
    config['global_settings']['hg19_genome_fasta'] = '/somewhere/else/genome.fasta'
    fs.create_symlink('/somewhere/else/genome.fasta', '/relative/genome.fasta')
    assert expected == helpers.make_apptainer_args(config, None)

    # Test default cache paths for array definition files
    fs.create_dir(cache_path)
    config['array_definition']['ExampleArray']['penncnv_pfb_file'] = '__cache-default__'
    # file not existing, binds write path
    expected = get_expected(expected_extra + [
        f"'{cache_path}/array_definitions/ExampleArray':'/outside/writearray/ExampleArray'"
    ])
    assert expected == helpers.make_apptainer_args(config, cache_path, not_existing_ok=True)
    # existing file directly bound
    fs.create_file(f'{cache_path}/array_definitions/ExampleArray/PennCNV-PFB_hg19.pfb')
    expected = get_expected(expected_extra + [
        f"'{cache_path}/array_definitions/ExampleArray/PennCNV-PFB_hg19.pfb':'/outside/ExampleArray/PennCNV-PFB_hg19.pfb'"
    ])
    assert expected == helpers.make_apptainer_args(config, cache_path)

    # test with cache_path & auto-creation of global paths
    config['global_settings'] = {
        'hg19_genome_fasta': '__default-ensemble__',
        'hg19_gtf_file': '__default-gencode__',
        'hg19_genomeInfo_file': '__default-UCSC__',
        'hg19_mehari_transcript_db': '__cache-default__',
        'dosage_sensitivity_scores': '__cache-default__',
    }
    expected_extra = []
    for global_file in ('fasta', 'gtf', 'genome_info', 'mehari_txdb', 'dosage_scores'):
        static_file = helpers.get_global_file(global_file, 'hg19', config['global_settings'], cache_path)
        fs.create_file(static_file)
        expected_extra += [f"'{static_file}':'/outside/static/{os.path.basename(static_file)}'"]
    expected_extra += [
        "'/relative_array/bpm_manifest.bpm':'/outside/ExampleArray/bpm_manifest.bpm'",
        f"'{cache_path}/array_definitions/ExampleArray/PennCNV-PFB_hg19.pfb':'/outside/ExampleArray/PennCNV-PFB_hg19.pfb'"
    ]
    assert get_expected(expected_extra) == helpers.make_apptainer_args(config, cache_path)

    # test with multiple arrays/genome versions
    config['global_settings'] .update({
        'hg38_genome_fasta': '__default-ensemble__',
        'hg38_gtf_file': '__default-gencode__',
        'hg38_genomeInfo_file': '__default-UCSC__',
        'hg38_mehari_transcript_db': '__cache-default__',
    })
    config['array_definition'].update({
        'ExampleArray2': {'genome_version': 'hg38'}
    })
    expected_extra = []
    for global_file in ('fasta', 'gtf', 'genome_info', 'mehari_txdb', 'dosage_scores'):
        static_file = helpers.get_global_file(global_file, 'hg38', config['global_settings'], cache_path)
        if not os.path.exists(static_file):
            fs.create_file(static_file)
        expected_extra += [f"'{static_file}':'/outside/static/{os.path.basename(static_file)}'"]
        static_file = helpers.get_global_file(global_file, 'hg19', config['global_settings'], cache_path)
        expected_extra += [f"'{static_file}':'/outside/static/{os.path.basename(static_file)}'"]
    expected_extra += [
        "'/relative_array/bpm_manifest.bpm':'/outside/ExampleArray/bpm_manifest.bpm'",
        f"'{cache_path}/array_definitions/ExampleArray/PennCNV-PFB_hg19.pfb':'/outside/ExampleArray/PennCNV-PFB_hg19.pfb'"
    ]
    assert get_expected(expected_extra) == helpers.make_apptainer_args(config, cache_path)

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

@pytest.fixture
def array_config():
    return textwrap.dedent(
        """\
        array_definition:
            ExampleArray:
                key: value
                key_file: filepath
        """
    )


def test_load_config(user_config, default_config, array_config, fs):

    fs.create_file('config.yaml', contents=user_config)
    fs.create_file(
        importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml'),
        contents=default_config
    )

    fake_args = MagicMock(config='config.yaml', no_cache=False, cache_path='/path/to/cache')

    config_dict = {
        'data_path': '/path/to/data',
        'cat': {'key': 'value',
                'nested': {'key': 'value2'}}}
    config_from_defaults = {
        'data_path': '/path/to/data',
        'cat': {'key': 'value',
                'nested': {'key': 'value2',
                           'key2': 'val3'}}}

    assert config_dict == helpers.load_config(fake_args, inbuilt_defaults=False)
    assert config_from_defaults == helpers.load_config(fake_args)

    # Test loading of array definition from present global config
    fs.create_file('/path/to/cache/global_array_definitions.yaml', contents=array_config)
    array_block = {'array_definition': {'ExampleArray': {'key': 'value', 'key_file': 'filepath'}}}
    config_from_defaults.update(array_block)
    assert config_from_defaults == helpers.load_config(fake_args)
    # # test overwrite of non-resolved file paths in user config with values from global config
    # fs.create_file(
    #     'config2.yaml',
    #     contents=user_config + '\n' + array_config.replace('filepath', 'non-existing-file')
    # )
    # fake_args.config = 'config2.yaml'
    # assert config_from_defaults == helpers.load_config(fake_args)

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
    sampletable = fs.create_file('sample_table.tsv', contents=sample_table_extra_cols)
    sample_data_df = helpers.read_sample_table('sample_table.tsv')

    all_ids = ['Cellline-A-MB', 'Cellline-A-WB', 'Cellline-B-MB', 'Cellline-B-1-cl1']

    # Test finding sample_ids based on matching entries in given column
    # The search sample ID itself is *never* included in the result, the reference sample always is
    assert set(all_ids[1:]) == helpers.collect_SNP_cluster_ids(
        'Cellline-A-MB',
        {'match_columns': ['Sample_Group'], 'id_columns': [], 'sample_ids': [], 'max_number_samples': 20},
        sample_data_df
    )
    assert {all_ids[0]} == helpers.collect_SNP_cluster_ids(
        'Cellline-A-WB',
        {'match_columns': ['Chip_Pos'], 'id_columns': [], 'sample_ids': [], 'max_number_samples': 20},
        sample_data_df)
    with pytest.raises(ConfigValueError):
        helpers.collect_SNP_cluster_ids(
            'Cellline-A-MB',
            {'match_columns': ['NonExisting'], 'id_columns': [], 'sample_ids': [], 'max_number_samples': 20},
            sample_data_df
        )
    # Find only reference sample
    assert {all_ids[0]} == helpers.collect_SNP_cluster_ids(
        'Cellline-A-WB',
        {'match_columns': [], 'id_columns': [], 'sample_ids': [], 'max_number_samples': 20},
        sample_data_df
    )
    # also needs work when no other samples are matching
    assert set() == helpers.collect_SNP_cluster_ids(
        'Cellline-A-MB',
        {'match_columns': ['Sample_Name'], 'id_columns': [], 'sample_ids': [], 'max_number_samples': 20},
        sample_data_df
    )
    # Do not use empty strings as matching value
    sample_data_df['empty_col'] = ''
    assert set() == helpers.collect_SNP_cluster_ids(
        'Cellline-A-MB',
        {'match_columns': ['empty_col'], 'id_columns': [], 'sample_ids': [], 'max_number_samples': 20},
        sample_data_df
    )
    # Test finding sample_ids based on values in given column
    assert set(all_ids[1:3]) == helpers.collect_SNP_cluster_ids(
        'Cellline-A-MB',
        {'match_columns': [], 'id_columns': ['Test_col'], 'sample_ids': [], 'max_number_samples': 20},
        sample_data_df
    )
    # reference & ID column
    assert {all_ids[0], all_ids[2]} == helpers.collect_SNP_cluster_ids(
        'Cellline-A-WB',
        {'match_columns': [], 'id_columns': ['Test_col'], 'sample_ids': [], 'max_number_samples': 20},
        sample_data_df
    )
    with pytest.raises(ConfigValueError):
        helpers.collect_SNP_cluster_ids(
            'Cellline-A-WB',
            {'match_columns': [], 'id_columns': ['Testcol'], 'sample_ids': [], 'max_number_samples': 20},
            sample_data_df
        )
    # test no/empty entry in given column
    assert set() == helpers.collect_SNP_cluster_ids(
        'Cellline-B-MB',
        {'match_columns': [], 'id_columns': ['Reference_Sample'], 'sample_ids': [], 'max_number_samples': 20},
        sample_data_df
    )
    # Test using sample_ids directly (+ ref sample)
    assert set([all_ids[0]] + all_ids[2:]) == helpers.collect_SNP_cluster_ids(
        'Cellline-A-WB',
        {'match_columns': [], 'id_columns': ['Test_col'], 'sample_ids': ['Cellline-B-1-cl1'], 'max_number_samples': 20},
        sample_data_df
    )
    # fail on any extracted entries that contain non_existing samples
    with pytest.raises(SampleConstraintError):
        helpers.collect_SNP_cluster_ids(
            'Cellline-B-MB',
            {'match_columns': [], 'id_columns': ['Test_col'], 'sample_ids': [], 'max_number_samples': 20},
            sample_data_df
        )
    #FIXME: add checks for log messages to the next 2 tests

    # Test exclusion of samples from a different array
    sampletable.set_contents(
        sample_table_extra_cols +
        "DifferentArray\tCellline-X-MB\t0\tMale\tCellline-X\t\tExampleCellines\t773456789000\tR01C03\tCellline-X-MB"
    )
    assert set(all_ids[1:]) == helpers.collect_SNP_cluster_ids(
        'Cellline-A-MB',
        {'match_columns': ['Sample_Group'], 'id_columns': [], 'sample_ids': [], 'max_number_samples': 20},
        sample_data_df
    )
    # test with max_number_samples
    assert {all_ids[3], all_ids[1]} == helpers.collect_SNP_cluster_ids(
        'Cellline-A-MB',
        {'match_columns': [], 'id_columns': ['Test_col'], 'sample_ids': ['Cellline-B-1-cl1'], 'max_number_samples': 2},
        sample_data_df
    )



def test_get_cache_dir(caplog, fs):
    config = {
        'global_settings': {
            'cache_dir': '/tmp/stemcnv_cache'
        }
    }
    args = MagicMock(
        no_cache=True,
        cache_path=None
    )
    fs.add_real_file(
        importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml'),
    )
    
    # Test with args.no_cache
    assert helpers.get_cache_dir(args, config) is None

    args.no_cache = False
    # test if set cache_path is not accessible
    args.cache_path = '/unwritable/cache'
    fs.create_dir('/unwritable/cache', perm_bits=444)
    assert helpers.get_cache_dir(args, config) is None
    assert [log.levelname for log in caplog.records[-3:]] == ['DEBUG', 'DEBUG', 'WARNING']
    assert "Failed to create or access cache in: /unwritable/cache" == caplog.records[-3].message
    assert "Could not use cache directory '/unwritable/cache'" == caplog.records[-2].message

    # Test with usable ste cache_path from args
    args.cache_path = '/writable/cache'
    assert helpers.get_cache_dir(args, config) == Path('/writable/cache')
    assert fs.get_object('/writable/cache') is not None

    # Test without set cache_path = fallback to config value
    args.cache_path = ''
    assert helpers.get_cache_dir(args, config) == Path('/tmp/stemcnv_cache')
    assert fs.get_object('/tmp/stemcnv_cache') is not None

    # Test fallback to default home dir path if config dir is unaccessible
    fs.remove_object('/tmp/stemcnv_cache')
    fs.create_dir('/tmp/stemcnv_cache', perm_bits=000)
    actual = helpers.get_cache_dir(args, config)
    assert actual == Path('~/.cache/stemcnv-check').expanduser()
    assert fs.get_object(Path('~/.cache/stemcnv-check').expanduser()) is not None
    assert [log.levelname for log in caplog.records[-3:]] == ['DEBUG', 'WARNING', 'DEBUG']
    assert "Cache path exists but is not accessible: /tmp/stemcnv_cache" == caplog.records[-3].message
    assert "Could not use cache directory '/tmp/stemcnv_cache'" == caplog.records[-2].message
    assert f"Created cache directory: {actual}" == caplog.records[-1].message

    # Test fallback top default cache path if config has no entry
    config = {}
    assert helpers.get_cache_dir(args, config) == Path('~/.cache/stemcnv-check').expanduser()


def test_get_global_file(fs):
    from stemcnv_check import ENSEMBL_RELEASE, MEHARI_DB_VERSION

    # Test with pre-defined files
    global_settings = OrderedDict({
        'hg19_genome_fasta': 'relative/genome.fasta',
        'hg19_gtf_file': 'relative/gtf.gtf',
        'hg19_genomeInfo_file': 'relative/genome.info',
        'hg19_mehari_transcript_db': 'relative/mehari_db',
    })
    for gtype, expected in zip(('fasta', 'gtf', 'genome_info', 'mehari_txdb'), global_settings.values()):
        assert helpers.get_global_file(gtype, 'hg19', global_settings, None) == expected

    # Test default cache paths
    global_settings = {
        'hg19_genome_fasta': '__default-ensemble__',
        'hg19_gtf_file': '__default-gencode__',
        'hg19_genomeInfo_file': '__default-UCSC__',
        'hg19_mehari_transcript_db': '__cache-default__',
        'hg38_genome_fasta': '__default-ensemble__',
        'hg38_gtf_file': '__default-gencode__',
        'hg38_genomeInfo_file': '__default-UCSC__',
        'hg38_mehari_transcript_db': '__cache-default__',
    }

    # Error without existing cache
    with pytest.raises(CacheUnavailableError):
        helpers.get_global_file('fasta', 'hg19', global_settings, None)

    # Test the predefined defaults with cache
    cache = '/path/to/cache'
    fs.create_dir(cache)

    expected = {
        'fasta': os.path.join(cache, 'fasta', 'homo_sapiens', f'{ENSEMBL_RELEASE}_{{genome}}', 'Homo_sapiens.{genome}.dna.primary_assembly.fa.gz'),
        'gtf': os.path.join(cache, 'static-data', 'gencode.{genome}.v45.gtf.gz'),
        'genome_info': os.path.join(cache, 'static-data', 'UCSC_{genome}_chromosome-info.tsv'),
        'mehari_txdb': os.path.join(cache, 'mehari-db', "mehari-data-txs-{genome}-ensembl-{MEHARI_DB_VERSION}.bin.zst")
    }
    # fasta and mehari need specific format versions!
    genome_format = defaultdict(
        lambda: {'hg38': 'hg38', 'hg19': 'hg19'},
        fasta={'hg38': 'GRCh38', 'hg19': 'GRCh37'},
        mehari_txdb={'hg38': 'GRCh38', 'hg19': 'GRCh37'}
    )

    for genome in ('hg38', 'hg19'):
        expected['fasta'] = os.path.join(cache, 'fasta', 'homo_sapiens', f'{ENSEMBL_RELEASE}_{{genome}}', 'Homo_sapiens.{genome}.dna.primary_assembly.fa.gz')
        for gtype in ('fasta', 'gtf', 'genome_info', 'mehari_txdb'):
            assert helpers.get_global_file(gtype, genome, global_settings, cache, fill_wildcards=False) == expected[gtype]
            assert helpers.get_global_file(gtype, genome, global_settings, cache) == expected[gtype].format(
                genome=genome_format[gtype][genome],
                MEHARI_DB_VERSION=MEHARI_DB_VERSION
            )


def test_get_array_file(fs):
    cache_path = '/path/to/cache'
    fs.create_dir(cache_path)
    config = {
        'global_settings': {
            'cache_dir': cache_path
        },
        'array_definition': {
            'ExampleArray': {
                'genome_version': 'hg19',
                'bpm_manifest': 'relative_array/bpm_manifest.bpm',
                'penncnv_pfb_file': 'relative_array/file.pfb',
                'penncnv_GCmodel_file': '__cache-default__',
                'array_density_file': '__cache-default__',
                'array_gaps_file': '__cache-default__',
            }
        }
    }
    expected = {
        'penncnv_pfb_file': 'relative_array/file.pfb',
        'penncnv_GCmodel_file': '/path/to/cache/array_definitions/ExampleArray/PennCNV-GCmodel_hg19.gcmodel',
        'array_density_file': '/path/to/cache/array_definitions/ExampleArray/density_hg19.bed',
        'array_gaps_file': '/path/to/cache/array_definitions/ExampleArray/gaps_hg19.bed',
    }
    # test expected paths
    for filekey, expected_file in expected.items():
        assert expected_file == helpers.get_array_file(filekey, 'ExampleArray', config, cache_path)

    # Test Cache error if no cache given, raised ONLY if file should be in cache
    assert expected['penncnv_pfb_file'] == helpers.get_array_file(
        'penncnv_pfb_file', 'ExampleArray', config, None
    )
    with pytest.raises(CacheUnavailableError):
        helpers.get_array_file('penncnv_GCmodel_file', 'ExampleArray', config, None)

    # Test ValueError if file not one of the array/penncnv files
    with pytest.raises(ValueError):
        helpers.get_array_file('bpm_manifest', 'ExampleArray', config, cache_path)
