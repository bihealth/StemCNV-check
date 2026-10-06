import importlib.resources
import ruamel.yaml as ruamel_yaml
from argparse import Namespace
from unittest.mock import patch
from stemcnv_check.__main__ import main
from stemcnv_check.app.check_input import check_config
from stemcnv_check import STEM_CNV_CHECK

@patch('stemcnv_check.__main__.logging')
def test_setup_files_config(mock_logging, sample_table_minimal, fs):

    mock_logging.add.return_value = 0
    mock_logging.remove.return_value = 0

    fs.add_real_file(importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'label_name_definitions.yaml'))
    fs.add_real_file(importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml'))
    fs.add_real_file(importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'sample_table_example.tsv'))
    fs.add_real_file(importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'allowedvalues_config.yaml'))
    fs.create_file('bpm_manifest.bpm')
    fs.create_file('pfb_file.pfb')
    fs.create_file('egt_file.egt')
    yaml = ruamel_yaml.YAML(typ='safe')
    #test_config = fs.create_file('test_config.yaml')

    for detail_level in ('minimal', 'medium', 'advanced', 'complete'):
        # Create config at different levels (+ sample_table)
        main(['setup-files', '--overwrite', '--config-details', detail_level])

        # Ensure config can be parsed
        with open('config.yaml') as f:
            config = yaml.load(f)
        # and validated (with added required user entries)!
        config['array_definition'] = {
            'ExampleArray': {
                'genome_version': 'hg19',
                'bpm_manifest_file': 'bpm_manifest.bpm',
                'egt_cluster_file': 'egt_file.egt',
                'csv_manifest_file': '',
                'penncnv_pfb_file': '__cache-default__',
                'penncnv_GCmodel_file': '__cache-default__',
                'array_density_file': '__cache-default__',
                'array_gaps_file': '__cache-default__',
            }
        }
        config['raw_data_folder'] = 'dummy'

        with open('test_config.yaml', 'w') as f:
            yaml.dump(config, f)
        check_config(
            Namespace(
                config='test_config.yaml',
                sample_table='sample_table.tsv',
                column_remove_regex=None,
                no_cache=False,
                cache_path='/cache'
            )
        )


    # with pytest.raises(FileNotFoundError) as e:
    #     ret_val = main(['run'])
    #     assert ret_val != 0
    #
    # with pytest.raises(SystemExit) as e:
    #     main(["--help"])
    #     assert e.value.code == 0
    #
    # fs.create_file(
    #     'sample_table.tsv',
    #     contents='Array_Name\tSample_ID\tChip_Name\tChip_Pos\tSex\tReference_Sample\tRegions_of_Interest\tSample_Group\n'
    #              'default\tSample1\t123456789000\tR01C03\tM\t\t\t\n'
    # )
    # fs.create_file('config.yaml', contents='array_definition: {"default" : {}}\n')
    # fs.create_dir('rundirectory')
    # default_config = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
    # allowed_values = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'allowedvalues_config.yaml')
    # fs.add_real_file(default_config, read_only=True)
    # fs.add_real_file(allowed_values, read_only=True)
    #
    # ret_val = main(['run'])
    # assert ret_val == 0
