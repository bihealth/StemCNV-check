import pytest
import textwrap
from loguru import logger
from _pytest.logging import LogCaptureFixture

from stemcnv_check.__main__ import filter_logs

@pytest.fixture
def caplog(caplog: LogCaptureFixture):
    logger.remove()
    handler_id = logger.add(
        caplog.handler,
        format="{message}",
        level=0,
        filter=lambda record: filter_logs(record) and record["level"].no >= caplog.handler.level and "MainThread" not in record["message"],
        enqueue=False,  # Set to 'True' if your test is spawning child processes.
        diagnose=False,
        backtrace=False
    )
    yield caplog
    logger.remove(handler_id)


@pytest.fixture
def sample_table_minimal():
    return textwrap.dedent(
        """\
        Array_Name\tChip_Name\tChip_Pos\tSample_ID\tSex\tReference_Sample\tRegions_of_Interest\tSample_Group
        ExampleArray\t123456789000\tR01C01\tCellline-A-MB\tFemale\t\tExample1|chr1:100000-200000\tExampleCellines
        ExampleArray\t123456789001\tR01C01\tCellline-A-WB\tFemale\tCellline-A-MB\tExample1|chr1:100000-200000;chr11:60000-70000\tExampleCellines
        """
    )