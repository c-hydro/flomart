"""
Library Features:

Name:          lib_utils_generic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220208'
Version:       '1.0.0'
"""

#######################################################################################
# Library
import logging
import numpy as np

from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
#######################################################################################


# -------------------------------------------------------------------------------------
# Method to pad or truncate list
def pad_or_truncate_list(some_list, target_len):
    return some_list[:target_len] + [0]*(target_len - len(some_list))
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to reduce dictionary 2 string
def reduce_dict_2_lists(data_dict):
    list_key, list_value = [], []
    for data_key, data_value in data_dict.items():
        list_key.append(data_key)
        list_value.append(data_value)

    return list_key, list_value
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to get nested value
def get_dict_nested_value(input_dict, nested_key):
    internal_dict_value = input_dict
    for k in nested_key:
        internal_dict_value = internal_dict_value.get(k, None)
        if internal_dict_value is None:
            return None
    return internal_dict_value
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to get recursively dictionary value
def get_dict_value(d, key, value=[]):

    for k, v in iter(d.items()):
        if isinstance(v, dict):
            if k == key:
                for kk, vv in iter(v.items()):
                    temp = [kk, vv]
                    value.append(temp)
            else:
                vf = get_dict_value(v, key, value)
                if isinstance(vf, list):
                    if vf:
                        vf_end = vf[0]
                    else:
                        vf_end = None
                elif isinstance(vf, np.ndarray):
                    vf_end = vf.tolist()
                else:
                    vf_end = vf
                if (vf_end is not None) and (vf_end not in value):

                    if not isinstance(vf_end, bool):
                        if vf_end:
                            if isinstance(value, list):
                                value.append(vf_end)
                            elif isinstance(value, str):
                                value = [value, vf_end]
                        else:
                            pass
                    else:
                        value.append(vf_end)
                else:
                    pass
        else:
            if k == key:
                if isinstance(v, np.ndarray):
                    value = v
                else:
                    value = v
    return value
# -------------------------------------------------------------------------------------
