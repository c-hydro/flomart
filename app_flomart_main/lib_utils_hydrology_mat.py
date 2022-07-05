"""
Library Features:

Name:          lib_utils_hydrology_mat
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220208'
Version:       '1.0.0'
"""

#######################################################################################
# Libraries
import logging
import os
import re

from datetime import datetime
from scipy.io import loadmat

import numpy as np
import pandas as pd

from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)
#######################################################################################


# -------------------------------------------------------------------------------------
# Create file hydro tag
def create_file_hydro_tag(section_ts_start, section_ts_end, section_ens=None, section_name=None,
                          time_format='%Y%m%d%H%M', tag_sep=':'):

    if (section_ens is not None) and (section_name is None):
        section_tag = section_ts_start.strftime(time_format) + '_' + section_ts_end.strftime(time_format) + \
                      tag_sep + section_ens
    elif (section_name is not None) and (section_ens is None):
        section_tag = section_ts_start.strftime(time_format) + '_' + section_ts_end.strftime(time_format) + \
                      tag_sep + section_name
    elif (section_name is None) and (section_ens is None):
        section_tag = section_ts_start.strftime(time_format) + '_' + section_ts_end.strftime(time_format)
    else:
        log_stream.error(' ===> File tag is not correctly defined')
        raise NotImplementedError('Case not implemented yet')

    return section_tag
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Parse file hydro name
def parse_file_hydro_name(file_name, file_data='observed'):

    file_parts = re.findall(r'\d+', file_name)
    file_root, file_extension = os.path.splitext(file_name)

    if not file_extension.endswith('mat'):
        log_stream.error(' ===> File extension must be "mat"')
        raise IOError('File type must be in ascii format')

    if file_parts.__len__() == 3:

        if file_data == 'observed':
            file_part_datetime_start = datetime.strptime(file_parts[0], "%Y%m%d%H%M")
            file_part_datetime_end = datetime.strptime(file_parts[1], "%Y%m%d%H%M")
            file_part_mask = None
            file_part_n_ens = None
        else:
            log_stream.error(
                ' ===> Parser of filename ' + file_name + ' failed for unknown type of file parts equal to 3')
            raise NotImplementedError('Case not implemented yet')
    else:
        log_stream.error(' ===> Parser of filename ' + file_name + ' failed for unknown format')
        raise NotImplementedError('Case not implemented yet')

    file_part_timestamp_start = pd.Timestamp(file_part_datetime_start)
    file_part_timestamp_end = pd.Timestamp(file_part_datetime_end)

    return file_part_timestamp_start, file_part_timestamp_end, file_part_mask, file_part_n_ens
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read hydro observed file in mat and txt format
def read_file_hydro_obs(section_name, file_name,
                        column_time_in='a1sDateVet', column_discharge_in='a1dQOss', column_level_in='a1dLivOssMean',
                        column_time_out='time', column_discharge_out='discharge', column_level_out='water_level'):

    file_data = loadmat(file_name)

    if column_time_in in list(file_data.keys()):
        file_time = file_data[column_time_in]
    else:
        log_stream.error(' ===> File column "' + column_time_in + '" observed not available in the datasets')
        raise IOError('Check your input file "' + file_name + '" to control the available fields')

    if column_discharge_in in list(file_data.keys()):
        file_discharge = file_data[column_discharge_in]
    else:
        log_stream.error(' ===> File column "' + column_discharge_in + '" observed not available in the datasets')
        raise IOError('Check your input file "' + file_name + '" to control the available fields')

    if column_level_in is not None:
        if column_level_in in list(file_data.keys()):
            file_water_level = file_data[column_level_in]
        else:
            log_stream.warning(' ===> File column "' + column_level_in + '" observed not available in the datasets')
            file_water_level = np.zeros(shape=(file_discharge.shape[0], 1))
            file_water_level[:, :] = -9999.0
    else:
        log_stream.warning(' ===> File column for water level observed variable is undefined')
        file_water_level = np.zeros(shape=(file_discharge.shape[0], 1))
        file_water_level[:, :] = -9999.0

    time_list = file_time[:, 0].tolist()
    discharge_list = file_discharge[:, 0].tolist()
    water_level_list = file_water_level[:, 0].tolist()

    section_period = []
    section_data_discharge = []
    section_data_water_level = []
    for time_step, discharge_step, water_level_step in zip(time_list, discharge_list, water_level_list):

        section_time = pd.Timestamp(pd.Timestamp(str(time_step[0])))
        section_point_discharge = float(discharge_step)
        section_point_water_level = float(water_level_step)

        section_period.append(section_time)
        section_data_discharge.append(section_point_discharge)
        section_data_water_level.append(section_point_water_level)

    section_series_discharge = pd.Series(index=section_period, data=section_data_discharge)
    section_series_water_level = pd.Series(index=section_period, data=section_data_water_level)

    return section_series_discharge, section_series_water_level

# -------------------------------------------------------------------------------------
