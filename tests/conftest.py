import pytest
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