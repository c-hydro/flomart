"""
Library Features:

Name:          lib_utils_hydrology_json
Author(s):     Matteo Darienzo and Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220711'
Version:       '2.0.0'
"""

#######################################################################################
# Libraries
import logging
import os
import re
from copy import deepcopy
from datetime import datetime
import numpy as np
import pandas as pd
import json
from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)


#######################################################################################


# -------------------------------------------------------------------------------------
# Method to read hydro simulated/observed file in json format
def read_file_hydro(section_name, file_name, variables_names):
    with open(file_name) as json_file:
        file_data = json.load(json_file)

    file_keys = list(file_data.keys())

    var_time = variables_names['time']
    var_discharge = variables_names['discharge']
    var_wlevel = variables_names['water_level']

    if var_time in list(file_keys):
        file_time = file_data[var_time]
        file_time = list(file_time.split(","))  # convert string to list of strings
    else:
        log_stream.error(' ===> File column "' + var_time + '" simulated not available in the datasets')
        raise IOError('Check your input file "' + file_name + '" to control the available fields')

    if var_discharge in list(file_keys):
        file_discharge = file_data[var_discharge]
        file_discharge = list(file_discharge.split(","))  # convert string to list of strings
    else:
        log_stream.error(' ===> File column "' + var_discharge + '" simulated not available in the datasets')
        raise IOError('Check your input file "' + file_name + '" to control the available fields')

    section_period = []
    section_data_discharge = []

    for time_step, discharge_step in zip(file_time, file_discharge):
        dt_tmp = pd.Timestamp(time_step)
        section_period.append(dt_tmp)

        section_point_discharge = float(discharge_step)
        section_data_discharge.append(section_point_discharge)

    section_series_discharge = pd.Series(index=section_period, data=section_data_discharge)

    # water level is not provided by json files!
    section_series_water_level = None

    return section_series_discharge, section_series_water_level


# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Create file hydro tag
def create_file_hydro_tag(section_ts_start, section_ts_end, section_ens=None, section_name=None,
                          time_format='%Y%m%d%H%M', tag_sep=':'):

    if (section_ens is not None) and (section_name is None):
        section_tag = section_ts_start.strftime(time_format) + '_' + section_ts_end.strftime(time_format) + \
                      tag_sep + section_ens
    elif (section_name is not None) and (section_ens is None):
        if section_ts_start == section_ts_end:
            section_tag = section_ts_end.strftime(time_format) + tag_sep + section_name
        else:
            section_tag = section_ts_start.strftime(time_format) + '_' + section_ts_end.strftime(time_format) + \
                          tag_sep + section_name
    elif (section_name is None) and (section_ens is None):
        if section_ts_start == section_ts_end:
            section_tag = section_ts_end.strftime(time_format)
        else:
            section_tag = section_ts_start.strftime(time_format) + '_' + section_ts_end.strftime(time_format)
    else:
        log_stream.error(' ===> File tag is not correctly defined')
        raise NotImplementedError('Case not implemented yet')

    return section_tag


# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Parse file hydro name
def parse_file_hydro_name(file_name, file_data='simulated'):

    match_parts = re.search(r'\d{4}\d{2}\d{2}\d{2}\d{2}', file_name)
    file_parts = [match_parts.group()]

    file_root, file_extension = os.path.splitext(file_name)

    if not file_extension.endswith('json'):
        log_stream.error(' ===> File extension must be "json"')
        raise IOError('File type must be in json format')

    if file_parts.__len__() == 1:

        # deterministic simulated file(s)
        if file_data == 'simulated':
            file_part_datetime_start = datetime.strptime(file_parts[0], "%Y%m%d%H%M")
            file_part_datetime_end = datetime.strptime(file_parts[0], "%Y%m%d%H%M")
            file_part_mask = None
            file_part_n_ens = None
        elif file_data == 'observed':
            file_part_datetime_start = datetime.strptime(file_parts[0], "%Y%m%d%H%M")
            file_part_datetime_end = datetime.strptime(file_parts[0], "%Y%m%d%H%M")
            file_part_mask = None
            file_part_n_ens = None
        else:
            log_stream.error(' ===> Parser of filename ' + file_name + ' failed for unknown file type')
            raise NotImplementedError('Case not implemented yet')
    else:
        log_stream.error(' ===> Parser of filename ' + file_name + ' failed for unknown format')
        raise NotImplementedError('Case not implemented yet')

    file_part_timestamp_start = pd.Timestamp(file_part_datetime_start)
    file_part_timestamp_end = pd.Timestamp(file_part_datetime_end)

    return file_part_timestamp_start, file_part_timestamp_end, file_part_mask, file_part_n_ens
# -------------------------------------------------------------------------------------
