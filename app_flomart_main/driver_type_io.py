"""
Class Features

Name:          driver_type_io
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231114'
Version:       '1.1.0'
"""

######################################################################################
# Library
import logging
import os
import pandas as pd

from lib_utils_hydrology_generic import analyze_obj_hydro, create_obj_hydro

from lib_utils_hydrology_generic import fill_obj_hydro

from lib_utils_hydrology_ascii import read_file_hydro_sim_drift as read_file_hydro_sim_ascii_drift
from lib_utils_hydrology_ascii import read_file_hydro_sim_hmc as read_file_hydro_sim_ascii_hmc
from lib_utils_hydrology_ascii import read_file_hydro_obs as read_file_hydro_obs_ascii
from lib_utils_hydrology_ascii import parse_file_hydro_name as parse_file_hydro_name_ascii
from lib_utils_hydrology_ascii import get_file_hydro_time as get_file_hydro_time_ascii
from lib_utils_hydrology_ascii import create_file_hydro_tag as create_file_hydro_tag_ascii

from lib_utils_hydrology_mat import read_file_hydro_obs as read_file_hydro_obs_mat
from lib_utils_hydrology_mat import parse_file_hydro_name as parse_file_hydro_name_mat
from lib_utils_hydrology_mat import create_file_hydro_tag as create_file_hydro_tag_mat

from lib_utils_hydrology_json import read_file_hydro as read_file_hydro_json
from lib_utils_hydrology_json import parse_file_hydro_name as parse_file_hydro_name_json
from lib_utils_hydrology_json import create_file_hydro_tag as create_file_hydro_tag_json

from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)
# Debug
# import matplotlib.pylab as plt
######################################################################################


# -------------------------------------------------------------------------------------
# Class DriverType
class DriverType:

    # -------------------------------------------------------------------------------------
    # Initialize class
    def __init__(self, section_name, file_name, file_time, variables_names,
                 section_id=-1,
                 file_type='json', data_type='simulated',
                 method_data_filling='ffill', method_data_occurrence='all'):

        self.section_name = section_name
        self.section_id = section_id
        self.file_name = file_name
        self.file_time = file_time
        self.variables_names = variables_names
        self.file_type = file_type
        self.data_type = data_type
        self.method_data_filling = method_data_filling
        self.method_data_occurrence = method_data_occurrence
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to get data
    def get_data(self):

        if (self.file_type == 'ascii_arpal_sim') or (self.file_type == 'ascii_arpal_sim_drift'):
            section_ts = self.wrap_data_ascii(
                self.section_name, self.section_id, self.file_name, self.file_time, file_data=self.data_type,
                file_method_filling=self.method_data_filling, file_method_occurrence=self.method_data_occurrence,
                file_type='drift')
        elif self.file_type == 'ascii_arpal_sim_hmc':
            section_ts = self.wrap_data_ascii(
                self.section_name, self.section_id, self.file_name, self.file_time, file_data=self.data_type,
                file_method_filling=self.method_data_filling, file_method_occurrence=self.method_data_occurrence,
                file_type='hmc')
        elif self.file_type == 'ascii_arpal_obs':
            section_ts = self.wrap_data_ascii(
                self.section_name, self.section_id, self.file_name, self.file_time, file_data=self.data_type,
                file_method_filling=self.method_data_filling, file_method_occurrence=self.method_data_occurrence)
        elif self.file_type == 'matlab':
            section_ts = self.wrap_data_mat(
                self.section_name, self.file_name, self.file_time, file_data=self.data_type,
                file_method_filling=self.method_data_filling, file_method_occurrence=self.method_data_occurrence)
        elif self.file_type == 'json':
            section_ts = self.wrap_data_json(
                self.section_name, self.file_name, self.file_time, variables_names=self.variables_names,
                file_data=self.data_type, file_method_filling=self.method_data_filling, file_method_occurrence=self.method_data_occurrence)
        else:
            log_stream.error(' ===> FileType "' + self.file_type + '" is not supported')
            raise NotImplementedError('Case not implemented yet')

        if section_ts is not None:
            if not isinstance(section_ts, pd.DataFrame):
                log_stream.error(' ===> Section TimeSeries obj is expected in DataFrame format')
                raise RuntimeError('Check your Section TimeSeries object to define it using the DataFrame format')

        else:
            log_stream.warning(' ===> Section TimeSeries is defined by NoneType. All files are not available')

        return section_ts
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # # Method to wrap ascii data
    @staticmethod
    def wrap_data_ascii(section_name, section_id, file_path_list, file_time,
                        file_data='simulated',
                        file_method_occurrence='first', file_method_filling='ffill', file_type='drift'):

        # check file list format
        if not isinstance(file_path_list, list):
            file_path_list = [file_path_list]

        # iterate over file list
        section_ts_discharge, section_ts_water_level, section_dframe_discharge = None, None, None
        for file_path_step in file_path_list:

            # get section series obj
            if file_data == 'simulated':
                if file_type == 'drift':
                    section_ts_discharge = read_file_hydro_sim_ascii_drift(section_name, file_path_step)
                elif file_type == 'hmc':
                    section_ts_discharge = read_file_hydro_sim_ascii_hmc(section_name, section_id, file_path_step)
                else:
                    log_stream.error(' ===> File type "' + file_type + '" is not supported')
                    raise NotImplementedError('Case not implemented yet')
            elif file_data == 'observed':
                section_ts_discharge, section_ts_water_level = read_file_hydro_obs_ascii(section_name, file_path_step)
            else:
                log_stream.error(' ===> File data "' + file_data + '" is not supported')
                raise NotImplementedError('Case not implemented yet')

            # check section series object
            if section_ts_discharge is not None:

                # get folder and file name(s)
                folder_name_step, file_name_step = os.path.split(file_path_step)
                # initialize dataframe object
                if section_dframe_discharge is None:
                    section_dframe_discharge = pd.DataFrame(index=file_time)

                # method to fill empty step(s)
                section_ts_discharge = fill_obj_hydro(
                    obj_hydro=section_ts_discharge, obj_method_filling=file_method_filling)
                # method to get file time(s)
                file_time_start, file_time_end = get_file_hydro_time_ascii(section_ts_discharge)

                # method to parse file name(s)
                file_part_timestamp_start, file_part_timestamp_end, \
                    file_part_mask, file_part_n_ens = parse_file_hydro_name_ascii(
                        file_name_step, file_data=file_data,
                        file_time_start=file_time_start, file_time_end=file_time_end)
                # method to create file tag(s)
                file_tag = create_file_hydro_tag_ascii(
                    file_part_timestamp_start, file_part_timestamp_end, file_part_n_ens, section_name)
                # organize file attributes
                section_ts_discharge.attrs = {
                    'file_name': file_name_step, 'section_name': section_name,
                    'file_mask': file_part_mask, 'file_ensemble': file_part_n_ens, 'file_tag': file_tag,
                    'time_start': file_part_timestamp_start, 'time_end': file_part_timestamp_end}
                # store data according to file tag
                section_dframe_discharge[file_tag] = section_ts_discharge

            # keep the occurrence of the first file or all file(s)
            if section_dframe_discharge is not None:
                if file_method_occurrence == 'first':
                    break
                elif file_method_occurrence == 'all':
                    pass
                else:
                    log_stream.error(' ===> File occurrence "' + file_method_occurrence + '" is not supported')
                    raise NotImplementedError('Case not implemented yet')

        # check the section collections object
        if section_dframe_discharge is None:
            log_stream.warning(' ===> Dataframe for section "' + section_name + '" is empty')

        return section_dframe_discharge
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to wrap mat data
    @staticmethod
    def wrap_data_mat(section_name, file_path_list, file_time,
                      file_data='observed',
                      file_method_occurrence='first', file_method_filling='interpolate'):

        if not isinstance(file_path_list, list):
            file_path_list = [file_path_list]

        section_ts_discharge, section_ts_water_level, section_dframe_discharge = None, None, None
        for file_path_step in file_path_list:

            if file_data == 'observed':
                section_ts_discharge, section_ts_water_level = read_file_hydro_obs_mat(section_name, file_path_step)
            else:
                log_stream.error(' ===> File data "' + file_data + '" is not supported')
                raise NotImplementedError('Case not implemented yet')

            if section_ts_discharge is not None:

                folder_name_step, file_name_step = os.path.split(file_path_step)

                if section_dframe_discharge is None:
                    section_dframe_discharge = pd.DataFrame(index=file_time)

                section_ts_discharge = fill_obj_hydro(
                    obj_hydro=section_ts_discharge, obj_method_filling=file_method_filling)

                file_part_timestamp_start, file_part_timestamp_end, \
                    file_part_mask, file_part_n_ens = parse_file_hydro_name_mat(
                        file_name_step, file_data=file_data)

                file_tag = create_file_hydro_tag_mat(
                    file_part_timestamp_start, file_part_timestamp_end, file_part_n_ens, section_name)

                section_ts_discharge.attrs = {
                    'file_name': file_name_step, 'section_name': section_name,
                    'file_mask': file_part_mask, 'file_ensemble': file_part_n_ens, 'file_tag': file_tag,
                    'time_start': file_part_timestamp_start, 'time_end': file_part_timestamp_end}

                section_dframe_discharge[file_tag] = section_ts_discharge

            if section_dframe_discharge is not None:
                if file_method_occurrence == 'first':
                    break
                elif file_method_occurrence == 'all':
                    pass
                else:
                    log_stream.error(' ===> File occurrence "' + file_method_occurrence + '" is not supported')
                    raise NotImplementedError('Case not implemented yet')

        if section_dframe_discharge is not None:
            section_dframe_discharge = fill_obj_hydro(
                obj_hydro=section_dframe_discharge, obj_method_filling=file_method_filling)
        else:
            log_stream.warning(' ===> Section "' + section_name + '" has not valid time series ')

        return section_dframe_discharge
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to wrap json data
    @staticmethod
    def wrap_data_json(section_name, file_path_list, file_time,
                       variables_names, file_data='simulated', file_method_occurrence='first',
                       file_method_filling='ffill'):

        if not isinstance(file_path_list, list):
            file_path_list = [file_path_list]

        section_ts_discharge, section_ts_water_level, section_dframe_discharge = None, None, None
        for file_path_step in file_path_list:

            if file_data == 'simulated':
                section_ts_discharge, section_ts_water_level = read_file_hydro_json(
                    section_name, file_path_step, variables_names)
            elif file_data == 'observed':
                section_ts_discharge, section_ts_water_level = read_file_hydro_json(
                    section_name, file_path_step, variables_names)
            else:
                log_stream.error(' ===> File data "' + file_data + '" is not supported')
                raise NotImplementedError('Case not implemented yet')

            if section_ts_discharge is not None:

                folder_name_step, file_name_step = os.path.split(file_path_step)

                if section_dframe_discharge is None:
                    section_dframe_discharge = pd.DataFrame(index=file_time)

                section_ts_discharge = fill_obj_hydro(
                    obj_hydro=section_ts_discharge, obj_method_filling=file_method_filling)

                file_part_timestamp_start, file_part_timestamp_end, \
                file_part_mask, file_part_n_ens = parse_file_hydro_name_json(
                     file_name_step, file_data=file_data)

                file_tag = create_file_hydro_tag_json(
                     file_part_timestamp_start, file_part_timestamp_end, file_part_n_ens, section_name)

                section_ts_discharge.attrs = {
                     'file_name': file_name_step, 'section_name': section_name,
                     'file_mask': file_part_mask, 'file_ensemble': file_part_n_ens, 'file_tag': file_tag,
                     'time_start': file_part_timestamp_start, 'time_end': file_part_timestamp_end}

                section_dframe_discharge[file_tag] = section_ts_discharge

            if section_dframe_discharge is not None:
                if file_method_occurrence == 'first':
                    break
                elif file_method_occurrence == 'all':
                    pass
                else:
                    log_stream.error(' ===> File occurrence "' + file_method_occurrence + '" is not supported')
                    raise NotImplementedError('Case not implemented yet')

        if section_dframe_discharge is not None:
            section_dframe_discharge = fill_obj_hydro(
                obj_hydro=section_dframe_discharge, obj_method_filling=file_method_filling)
        else:
            log_stream.warning(' ===> Section "' + section_name + '" has not valid time series ')

        return section_dframe_discharge
        # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
