import pytest
import importlib.resources
from unittest.mock import patch
from stemcnv_check.__main__ import main
from stemcnv_check import STEM_CNV_CHECK

@patch('stemcnv_check.__main__.check_config')
@patch('stemcnv_check.__main__.run_stemcnv_check_workflow')
def test_main_run(mock_check, mock_run_workflow, fs):

    mock_check.return_value = 0
    mock_run_workflow.return_value = 0

    with pytest.raises(FileNotFoundError) as e:
        main([])
        assert e.value.code != 0

    with pytest.raises(SystemExit) as e:
        main(["--help"])
        assert e.value.code == 0

    fs.create_file('sample_table.txt',
                   contents='Sample_ID\tChip_Name\tChip_Pos\tSex\tReference_Sample\nSample1\t123456789000\tR01C03\tM\t\n')
    fs.create_file('config.yaml', contents='static_data: {}\n')
    fs.create_dir('rundirectory')
    default_config = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
    allowed_values = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'allowedvalues_config.yaml')
    fs.add_real_file(default_config, read_only=True)
    fs.add_real_file(allowed_values, read_only=True)

    with pytest.raises(SystemExit) as e:
        main([])
        assert e.value.code == 0

