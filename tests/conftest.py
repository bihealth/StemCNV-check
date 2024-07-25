import pytest
from loguru import logger
from _pytest.logging import LogCaptureFixture

@pytest.fixture
def caplog(caplog: LogCaptureFixture):
    logger.remove()
    handler_id = logger.add(
        caplog.handler,
        format="{message}",
        level=0,
        filter=lambda record: record["level"].no >= caplog.handler.level and "MainThread" not in record["message"],
        enqueue=False,  # Set to 'True' if your test is spawning child processes.
        diagnose=False,
        backtrace=False
    )
    yield caplog
    logger.remove(handler_id)