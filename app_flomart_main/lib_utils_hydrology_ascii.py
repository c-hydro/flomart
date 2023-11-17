"""
Library Features:

Name:          lib_utils_hydrology_ascii
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231114'
Version:       '1.1.0'
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

from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)
#######################################################################################


# -------------------------------------------------------------------------------------
# method to read hydrological section file
def read_file_hydro_section(file_name, file_sep=' ', file_columns=None):

    if file_columns is None:
        file_columns = ['x', 'y',
                        'basin_name', 'section_name', 'org_name',
                        'drainage_area', 'thr_1', 'thr_2', 'domain_name', 'code']

    file_table = pd.read_table(file_name, sep=file_sep, names=file_columns)

    file_dict = {}
    for row_id, row_fields in enumerate(file_table.iterrows()):
        row_dict = row_fields[1].to_dict()
        row_dict['section_loc'] = int(row_id)

        file_dict['section_{:}'.format(row_id)] = row_dict

    return file_dict

# -------------------------------------------------------------------------------------


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
# method to get time(s)
def get_file_hydro_time(section_series):
    if isinstance(section_series, pd.Series):
        section_time_start = section_series.index[0].to_datetime64()
        section_time_end = section_series.index[-1].to_datetime64()
    else:
        section_time_start, section_time_end = None, None
    return section_time_start, section_time_end
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Parse file hydro name
def parse_file_hydro_name(file_name, file_data='simulated', file_time_start=None, file_time_end=None):

    file_parts = re.findall(r'\d+', file_name)
    file_root, file_extension = os.path.splitext(file_name)

    if not file_extension.endswith('txt'):
        log_stream.error(' ===> File extension must be "txt"')
        raise IOError('File type must be in ascii format')

    if file_parts.__len__() == 3:

        # probabilistic simulated file(s)
        if file_data == 'simulated':
            file_part_datetime_start = datetime.strptime(file_parts[0], "%y%j%H%M")
            file_part_datetime_end = datetime.strptime(file_parts[1][:-2], "%y%j%H%M")
            file_part_mask = file_parts[1][-2:]
            file_part_n_ens = file_parts[2]
        elif file_data == 'observed':
            file_part_datetime_start = datetime.strptime(file_parts[0], "%Y%m%d%H%M")
            file_part_datetime_end = datetime.strptime(file_parts[1], "%Y%m%d%H%M")
            file_part_mask = None
            file_part_n_ens = None
        else:
            log_stream.error(
                ' ===> Parser of filename ' + file_name + ' failed for unknown type of file parts equal to 3')
            raise NotImplementedError('Case not implemented yet')

    elif file_parts.__len__() == 2:

        # deterministic simulated file(s)
        if file_data == 'simulated':
            file_part_datetime_start = datetime.strptime(file_parts[0], "%y%j%H%M")
            file_part_datetime_end = datetime.strptime(file_parts[1][:-2], "%y%j%H%M")
            file_part_mask = file_parts[1][-2:]
            file_part_n_ens = None
        elif file_data == 'observed':
            file_part_datetime_start = datetime.strptime(file_parts[0], "%Y%m%d%H%M")
            file_part_datetime_end = datetime.strptime(file_parts[1], "%Y%m%d%H%M")
            file_part_mask = None
            file_part_n_ens = None

        else:
            log_stream.error(' ===> Parser of filename ' + file_name + ' failed for unknown file type')
            raise NotImplementedError('Case not implemented yet')

    elif file_parts.__len__() == 1:

        # deterministic simulated file(s)
        if file_data == 'simulated':
            if file_time_start is not None:
                file_part_datetime_start = file_time_start
            else:
                file_part_datetime_start = datetime.strptime(file_parts[0], "%Y%m%d%H%M")
            if file_time_end is not None:
                file_part_datetime_end = file_time_end
            else:
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


# -------------------------------------------------------------------------------------
# Method to read hydro observed file in mat and txt format
def read_file_hydro_obs(section_name, file_name,
                        column_time_in='a1sDateVet', column_discharge_in='a1dQOssMean', column_level_in='a1dLivOssMean',
                        column_time_out='time', column_discharge_out='discharge', column_level_out='water_level'):

    file_table = pd.read_table(file_name, sep=' ')

    if column_time_in in list(file_table.columns):
        file_time = file_table[column_time_in].values
    else:
        log_stream.error(' ===> File column "' + column_time_in + '" observed not available in the datasets')
        raise IOError('Check your input file "' + file_name + '" to control the available fields')

    if column_discharge_in in list(file_table.columns):
        file_discharge = file_table[column_discharge_in].values
    else:
        log_stream.error(' ===> File column "' + column_discharge_in + '" observed not available in the datasets')
        raise IOError('Check your input file "' + file_name + '" to control the available fields')

    if column_level_in is not None:
        if column_level_in in list(file_table.columns):
            file_water_level = file_table[column_level_in].values
        else:
            log_stream.warning(' ===> File column "' + column_level_in + '" observed not available in the datasets')
            file_water_level = np.zeros(shape=file_discharge.shape[0])
            file_water_level[:] = -9999.0
    else:
        log_stream.warning(' ===> File column for water level observed variable is undefined')
        file_water_level = np.zeros(shape=file_discharge.shape[0])
        file_water_level[:] = -9999.0

    time_list = file_time.tolist()
    discharge_list = file_discharge.tolist()
    water_level_list = file_water_level.tolist()

    section_period = []
    section_data_discharge = []
    section_data_water_level = []
    for time_step, discharge_step, water_level_step in zip(time_list, discharge_list, water_level_list):

        section_time = pd.Timestamp(pd.Timestamp(str(time_step)))
        section_point_discharge = float(discharge_step)
        section_point_water_level = float(water_level_step)

        section_period.append(section_time)
        section_data_discharge.append(section_point_discharge)
        section_data_water_level.append(section_point_water_level)

    section_series_discharge = pd.Series(index=section_period, data=section_data_discharge)
    section_series_water_level = pd.Series(index=section_period, data=section_data_water_level)

    return section_series_discharge, section_series_water_level

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read hydro simulated file in ascii format (drift)
def read_file_hydro_sim_drift(section_name, file_name, column_time_idx=0):

    file_data = pd.read_table(file_name)
    file_data = file_data.loc[:, ~file_data.columns.str.match("Unnamed")]

    file_columns = list(file_data.columns)
    if isinstance(file_columns, list):
        if file_columns.__len__() == 1:
            file_cols_tmp = file_columns[0].split(' ')
            file_cols_filtered = list(filter(None, file_cols_tmp))
        elif file_columns.__len__() > 1:
            file_cols_filtered = deepcopy(file_columns)
        else:
            log_stream.error(' ===> File columns dimension is not allowed')
            raise NotImplementedError('Case not implemented yet')
    else:
        log_stream.error(' ===> File column format is not supported')
        raise NotImplementedError('Case not implemented yet')

    if section_name in file_cols_filtered:

        column_section_idx = file_cols_filtered.index(section_name)
        file_data_table = list(file_data.values)

        section_period = []
        section_data = []
        for file_data_row in file_data_table:

            file_row = list(file_data_row)
            if file_columns.__len__() == 1:
                file_data_parts = file_data_row[0].split(' ')
                file_data_parts = list(filter(None, file_data_parts))
            elif file_columns.__len__() > 1:
                file_data_parts = deepcopy(file_row)
            else:
                log_stream.error(' ===> File row dimension is not allowed')
                raise NotImplementedError('Case not implemented yet')

            section_time_string = str(int(file_data_parts[column_time_idx]))
            section_point = float(file_data_parts[column_section_idx])

            section_time_stamp = pd.Timestamp(section_time_string)

            section_period.append(section_time_stamp)
            section_data.append(section_point)

        section_series = pd.Series(index=section_period, data=section_data)

    else:
        section_series = None

    return section_series

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read hydro simulated file in ascii format (hmc)
def read_file_hydro_sim_hmc(section_name, section_id, file_name, column_time_idx=0):

    file_data = pd.read_table(file_name)
    file_data = file_data.loc[:, ~file_data.columns.str.match("Unnamed")]

    row_n = file_data.shape[0]

    time_obj, data_obj, name_obj = [], None, None
    for row_id, row_obj in enumerate(file_data.iterrows()):
        string_row = row_obj[1].values[0]
        parts_row_raw = string_row.split(' ')
        parts_row_def = [x for x in parts_row_raw if x]

        time_tmp = parts_row_def[0]
        time_stamp = pd.Timestamp(time_tmp)
        time_obj.append(time_stamp)

        row_tmp = parts_row_def[1:]
        row_data = [float(x) for x in row_tmp]
        row_data = np.asarray(row_data)
        col_n = row_data.shape[0]

        if data_obj is None:
            data_obj = np.zeros(shape=(row_n, col_n))
            data_obj[:, :] = np.nan
        if name_obj is None:
            name_obj = []
            for col_step in range(0, col_n):
                name_obj.append('section_{:}'.format(col_step))

        data_obj[row_id, :] = row_data

    data_df = pd.DataFrame(data=data_obj, index=time_obj)
    data_df.columns = name_obj
    data_df.index.name = 'time'

    tmp_series = data_df.iloc[:, section_id]
    section_period = pd.DatetimeIndex(tmp_series.index.values)
    section_data = tmp_series.values

    section_series = pd.Series(index=section_period, data=section_data)

    return section_series

# -------------------------------------------------------------------------------------
