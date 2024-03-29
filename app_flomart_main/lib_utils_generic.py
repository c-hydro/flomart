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
import warnings
import numpy as np
try:
    import h5py
except ImportError:
    warnings.warn(" ===> H5py library is not imported. File .mat will not correctly read")

from copy import deepcopy

from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
#######################################################################################


# -------------------------------------------------------------------------------------
# Method to convert array 2 list
def convert_array_2_list(data_in):
    data_out = None
    if not isinstance(data_in, list):
        if data_in.ndim == 1 or data_in.shape[0] == 1:
            for data in data_in[0]:

                if data_out is None:
                    data_out = []
                if isinstance(data, np.ndarray):
                    data_filter = data[0]

                    if isinstance(data_filter, np.ndarray):
                        if data_filter.shape[0] == 1:
                            data_filter = data_filter[0]
                        else:
                            log_stream.error(' ===> Error in getting data value')
                            raise NotImplementedError('Case not implemented yet')
                else:
                    data_filter = data

                if isinstance(data_filter, str):
                    data_tmp = str(data_filter)
                elif isinstance(data_filter, (int, np.integer)):
                    data_tmp = int(data_filter)
                elif isinstance(data_filter, (float, np.floating)):
                    data_tmp = float(data_filter)
                elif isinstance(data_filter, h5py._hl.dataset.Dataset):
                    print('ciao')
                else:
                    log_stream.error(' ===> Error in parsering data value')
                    raise NotImplementedError('Case not implemented yet')
                data_out.append(data_tmp)
        else:
            data_out = deepcopy(data_in)

    else:
        data_out = deepcopy(data_in)

    return data_out
# -------------------------------------------------------------------------------------


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
