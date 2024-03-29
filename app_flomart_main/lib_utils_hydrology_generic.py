"""
Library Features:

Name:          lib_utils_hydrology
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220208'
Version:       '1.0.0'
"""

#######################################################################################
# Libraries
import logging

from copy import deepcopy

import numpy as np
import pandas as pd

from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)
#######################################################################################


# -------------------------------------------------------------------------------------
# method to map description hydro and hydraulic
def map_description_hydro_2_hydraulic(hydro_data_in, hydraulic_data):

    if hydro_data_in is not None:
        section_data_out = {}
        for section_id, section_fields in hydro_data_in.items():
            section_description_tmp = section_fields['section_name']

            section_description_def, hydraulic_id = None, None
            for hydraulic_id, hydraulic_fields in hydraulic_data.items():

                hydraulic_description = hydraulic_fields['description']
                if 'alias' in list(hydraulic_fields.keys()):
                    hydraulic_alias = hydraulic_fields['alias']
                else:
                    hydraulic_alias = None

                if section_description_tmp == hydraulic_description:
                    section_description_def = deepcopy(section_description_tmp)
                else:
                    if hydraulic_alias is not None:
                        if section_description_tmp in hydraulic_alias:
                            idx_alias = hydraulic_alias.index(section_description_tmp)
                            tmp_alias = hydraulic_alias[idx_alias]
                            section_description_def = deepcopy(hydraulic_description)

                if section_description_def is not None:
                    break

            # fill the section fields
            if section_description_def is not None:
                # add fields to the hydro datasets
                section_fields['section_name'] = section_description_def
                # update the hydro tag according to the hydraulic tag
                section_data_out[hydraulic_id] = section_fields
            else:
                log_stream.warning(' ===> Section "' + section_description_tmp +
                                   '" is not available in the hydraulic description or alias')
    else:
        section_data_out = None

    if section_data_out is not None:
        section_data_out = dict(sorted(section_data_out.items()))

    return section_data_out
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to create obj hydro
def create_obj_hydro(file_index, file_list_key, file_list_data):
    file_data = {}
    for file_step_key, file_step_data in zip(file_list_key, file_list_data):
        file_data[file_step_key] = file_step_data
    file_dframe = pd.DataFrame(index=file_index, data=file_data)
    return file_dframe
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to analyze obj hydro
def analyze_obj_hydro(file_obj_in, file_source='simulated', file_method='max',
                      file_tag_discharge='discharge', file_tag_type='type'):
    file_obj_out = {}
    for file_key, file_dframe_in in file_obj_in.items():

        if file_dframe_in is not None:
            file_time = list(file_dframe_in.index)
            file_type = [file_source] * file_time.__len__()

            if file_method is not None:
                if file_method == 'max':
                    file_values = list(file_dframe_in.max(axis=1).values)
                elif file_method == 'mean':
                    file_values = list(file_dframe_in.mean(axis=1).values)
                else:
                    log_stream.error(' ===> Analyze method for file values is not correctly defined')
                    raise NotImplementedError('Case not implemented yet')
            else:

                if file_dframe_in.ndim == 1:
                    file_values = list(file_dframe_in.values)
                elif file_dframe_in.ndim > 1:
                    file_values = list(file_dframe_in.values[:, 0])
                else:
                    log_stream.error(' ===> Dataframe obj dimension(s) are not supported')
                    raise NotImplementedError('Case not implemented yet')

            file_dframe_out = create_obj_hydro(file_time,
                                               [file_tag_discharge, file_tag_type], [file_values, file_type])

        else:
            file_dframe_out = None

        file_obj_out[file_key] = file_dframe_out

    return file_obj_out

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to fill time series hydro
def fill_obj_hydro(obj_hydro, obj_method_filling='ffill', obj_limit_filling=10,
                   obj_value_min=0.0, obj_value_max=None, obj_nodata=-9999.0):

    if obj_value_min is not None:
        obj_hydro[obj_hydro < obj_value_min] = np.nan
    if obj_value_max is not None:
        obj_hydro[obj_hydro > obj_value_max] = np.nan

    time_valid_first = obj_hydro.first_valid_index()
    time_valid_last = obj_hydro.last_valid_index()

    if obj_method_filling is not None:
        if obj_method_filling == 'ffill':
            obj_hydro_filled = obj_hydro.ffill(axis='rows', limit=obj_limit_filling)
        elif obj_method_filling == 'bfill':
            obj_hydro_filled = obj_hydro.bfill(axis='rows', limit=obj_limit_filling)
        elif obj_method_filling == 'interpolate':
            obj_hydro_filled = obj_hydro.interpolate(method='values', limit=obj_limit_filling)
        else:
            log_stream.error(' ===> Method to fill data "' + obj_method_filling + '" is not suppoerted')
            raise NotImplementedError('Case not implemented yet')
    else:
        obj_hydro_filled = deepcopy(obj_hydro)

    obj_hydro_filled = obj_hydro_filled.fillna(obj_nodata)
    obj_hydro_filtered = obj_hydro_filled[time_valid_first:time_valid_last]

    return obj_hydro_filtered

# -------------------------------------------------------------------------------------
