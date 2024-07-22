# -*- coding: utf-8 -*-
"""Exception classes for pipeline"""


class SampleFormattingError(Exception):
    """Raised when the sample_table file contains misformatted values"""


class SampleConstraintError(Exception):
    """Raised when the sample_table file contains missing values or values not matching constraints"""


class ConfigValueError(Exception):
    """Raised when the config file contains values not allowed in a given field"""
