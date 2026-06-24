"""
Step 1: Extract sequence/structural features and BE-QA (hypothesis testing)

Usage example: python be3d_step1_hypothesis.py ./yaml/dnmt3a_local.yaml
"""

import yaml
import warnings

def get_required(config, key, expected_type):
    """
    Fetch a required parameter from config, raising an error if missing or wrong type. 
    """
    keys = key.split('.')
    value = config
    for k in keys:
        if not isinstance(value, dict) or k not in value:
            raise KeyError(f"Required parameter '{key}' not found in config.")
        value = value[k]
    if not isinstance(value, expected_type):
        type_label = (
            ' | '.join(t.__name__ for t in expected_type)
            if isinstance(expected_type, tuple)
            else expected_type.__name__
        )
        raise TypeError(
            f"Required parameter '{key}' expected {type_label}, "
            f"got {type(value).__name__} (value: {value!r})."
        )
    return value
 
def get_optional(config, key, expected_type, default):
    """
    Fetch an optional parameter from config, warning and defaulting if missing or wrong type.
    """
    keys = key.split('.')
    value = config
    for k in keys:
        if not isinstance(value, dict) or k not in value:
            warnings.warn(
                f"Optional parameter '{key}' not found in config. "
                f"Defaulting to {default!r}.",
                UserWarning, stacklevel=2,
            )
            return default
        value = value[k]
    if not isinstance(value, expected_type):
        type_label = (
            ' | '.join(t.__name__ for t in expected_type)
            if isinstance(expected_type, tuple)
            else expected_type.__name__
        )
        warnings.warn(
            f"Optional parameter '{key}' expected {type_label}, "
            f"got {type(value).__name__} (value: {value!r}). "
            f"Defaulting to {default!r}.",
            UserWarning, stacklevel=2,
        )
        return default
    return value

def load_config(config_yaml):
    with open(config_yaml, "r") as file:
        config = yaml.safe_load(file)
    return config

def find_union(input, pthr_str):
    if input[0] == f'p<{pthr_str}' or input[1] == f'p<{pthr_str}':
        return f'p<{pthr_str}'
    elif input[0] == f'p>={pthr_str}' or input[1] == f'p>={pthr_str}':
        return f'p>={pthr_str}'
    else:
        return '-'
        
