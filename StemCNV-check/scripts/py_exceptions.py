# -*- coding: utf-8 -*-
"""Exception classes for pipeline"""

class SampletableDefaultColError(Exception):
    """Raised when the sample_table file does not contain all necessary default columns"""

class SampleIDNonuniqueError(Exception):
    """Raised when the sample_table file contains non-unique Sample_ID values"""

class SampletableReferenceError(Exception):
    """Raised when the sample_table file contains a reference that does not exist or has a non-matching sex"""

class SampletableRegionError(Exception):
    """Raised when the sample_table file contains Regions_of_Interest entries that are not properly formatted"""

class SampleConstraintError(Exception):
    """Raised when the sample_table file contains values not matching the constraints (general or defined in config)"""

class InputFileError(Exception):
    """Raised when the static input files are not defined or missing"""

class ConfigValueError(Exception):
    """Raised when the config file contains values not allowed in a given field"""

class ConfigValueWarning(Warning):
    """Raised when the config file contains values not allowed in a given field"""

class ConfigKeyWarning(Warning):
    """Raised when the config file contains entries that are not defined or deprecated"""
