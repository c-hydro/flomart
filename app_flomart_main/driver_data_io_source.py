"""
Class Features

Name:          driver_data_io_source
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20200515'
Version:       '1.0.0'
"""

######################################################################################
# Library
import logging
import os
import numpy as np
import pandas as pd
import glob

from copy import deepcopy

from lib_utils_hydrology_generic import analyze_obj_hydro, create_obj_hydro
from lib_utils_io import read_obj, write_obj
from lib_utils_system import fill_tags2string, make_folder
from lib_utils_generic import get_dict_value, get_dict_nested_value
from lib_info_args import logger_name, time_format_algorithm

from driver_type_io import DriverType

# Logging
log_stream = logging.getLogger(logger_name)
# Debug
# import matplotlib.pylab as plt
######################################################################################


# -------------------------------------------------------------------------------------
# Class DriverDischarge
class DriverDischarge:

    # -------------------------------------------------------------------------------------
    # Initialize class
    def __init__(self, time_now, time_run, geo_data_collection, src_dict, anc_dict,
                 alg_ancillary=None, alg_template_tags=None,
                 flag_discharge_data_sim='discharge_data_simulated',
                 flag_discharge_data_obs='discharge_data_observed',
                 flag_cleaning_anc_discharge_sim=True, flag_cleaning_anc_discharge_obs=True):

        self.time_now = time_now
        self.time_run = time_run
        self.geo_data_collection = geo_data_collection

        self.flag_discharge_data_sim = flag_discharge_data_sim
        self.flag_discharge_data_obs = flag_discharge_data_obs

        self.alg_ancillary = alg_ancillary

        self.alg_template_tags = alg_template_tags
        self.file_name_tag = 'file_name'
        self.folder_name_tag = 'folder_name'
        self.type_tag = 'type'
        self.variables_tag = 'variables'
        self.method_data_occurrence_tag = 'method_data_occurrence'
        self.method_data_analysis_tag = 'method_data_analysis'
        self.method_data_filling_tag = 'method_data_filling'
        self.method_data_null_tag = 'method_data_null'
        self.time_period_tag = 'time_period'
        self.time_rounding_tag = 'time_rounding'
        self.time_frequency_tag = 'time_frequency'

        self.var_name_time = 'time'
        self.var_name_discharge = 'discharge'
        self.var_name_water_level = 'water_level'
        self.var_name_type = 'type'

        self.domain_name_list = self.alg_ancillary['domain_name']
        self.scenario_type = self.alg_ancillary['scenario_type']
        self.scenario_boundary = self.alg_ancillary['scenario_boundary']

        domain_section_dict = {}
        for domain_name_step in self.domain_name_list:
            domain_section_list = get_dict_value(geo_data_collection[domain_name_step], 'name_point_outlet', [])
            domain_section_dict[domain_name_step] = domain_section_list
        self.domain_section_dict = domain_section_dict

        domain_hydro_dict = {}
        for domain_name_step in self.domain_name_list:
            domain_hydro_list = get_dict_value(geo_data_collection[domain_name_step], 'name_point_obs', [])
            domain_hydro_dict[domain_name_step] = domain_hydro_list
        self.domain_hydro_dict = domain_hydro_dict

        domain_description_dict = {}  # MATTEO: add section description:
        for domain_name_step in self.domain_name_list:
            domain_description_list = get_dict_value(geo_data_collection[domain_name_step], 'section_description', [])
            domain_description_dict[domain_name_step] = domain_description_list
        self.domain_description_dict = domain_description_dict

        self.folder_name_discharge_sim = src_dict[self.flag_discharge_data_sim][self.folder_name_tag]
        self.file_name_discharge_sim = src_dict[self.flag_discharge_data_sim][self.file_name_tag]
        self.variables_discharge_sim = src_dict[self.flag_discharge_data_sim][self.variables_tag]
        self.type_sim = src_dict[self.flag_discharge_data_sim][self.type_tag]
        self.method_data_occurrence_sim = src_dict[self.flag_discharge_data_sim][self.method_data_occurrence_tag]
        self.method_data_analysis_sim = src_dict[self.flag_discharge_data_sim][self.method_data_analysis_tag]
        self.method_data_filling_sim = src_dict[self.flag_discharge_data_sim][self.method_data_filling_tag]
        self.method_data_null_sim = src_dict[self.flag_discharge_data_sim][self.method_data_null_tag]
        self.time_period_discharge_sim = src_dict[self.flag_discharge_data_sim][self.time_period_tag]
        self.time_rounding_discharge_sim = src_dict[self.flag_discharge_data_sim][self.time_rounding_tag]
        self.time_frequency_discharge_sim = src_dict[self.flag_discharge_data_sim][self.time_frequency_tag]

        self.folder_name_discharge_obs = src_dict[self.flag_discharge_data_obs][self.folder_name_tag]
        self.file_name_discharge_obs = src_dict[self.flag_discharge_data_obs][self.file_name_tag]
        self.variables_obs = src_dict[self.flag_discharge_data_obs][self.variables_tag]
        self.type_obs = src_dict[self.flag_discharge_data_obs][self.type_tag]
        self.method_data_occurrence_obs = src_dict[self.flag_discharge_data_obs][self.method_data_occurrence_tag]
        self.method_data_analysis_obs = src_dict[self.flag_discharge_data_obs][self.method_data_analysis_tag]
        self.method_data_filling_obs = src_dict[self.flag_discharge_data_obs][self.method_data_filling_tag]
        self.method_data_null_obs = src_dict[self.flag_discharge_data_obs][self.method_data_null_tag]
        self.time_period_discharge_obs = src_dict[self.flag_discharge_data_obs][self.time_period_tag]
        self.time_rounding_discharge_obs = src_dict[self.flag_discharge_data_obs][self.time_rounding_tag]
        self.time_frequency_discharge_obs = src_dict[self.flag_discharge_data_obs][self.time_frequency_tag]

        self.file_prefix_sim, self.file_sep_sim, self.file_elem_sim = None, None, None
        if (self.folder_name_discharge_sim is not None) and (self.file_name_discharge_sim is not None):
            self.file_path_discharge_sim = self.define_file_discharge(
                self.time_run, self.folder_name_discharge_sim, self.file_name_discharge_sim,
                file_name_prefix=self.file_prefix_sim, file_name_elem=self.file_elem_sim,
                file_name_sep=self.file_sep_sim,
                extra_args={'section_name_obj': self.domain_section_dict,
                            'section_description_obj': self.domain_description_dict     # MATTEO: add description
                            })
        else:
            log_stream.error(' ===> Source files of "overland_flow" are not defined ')
            raise IOError('Overflow datasets is needed by the application.')

        self.file_prefix_obs, self.file_sep_obs, self.file_elem_obs = None, None, None
        self.file_path_discharge_obs = self.define_file_discharge(
            self.time_run, self.folder_name_discharge_obs, self.file_name_discharge_obs,
            file_name_prefix=self.file_prefix_obs, file_name_elem=self.file_elem_obs, file_name_sep=self.file_sep_obs,
            extra_args={'section_name_obj': self.domain_hydro_dict,
                        'time_rounding': self.time_rounding_discharge_obs,
                        'time_frequency': self.time_frequency_discharge_obs,
                        'time_period': self.time_period_discharge_obs})

        self.var_time_discharge_sim, self.var_discharge_discharge_sim, \
            self.var_wlevel_discharge_sim = self.define_file_variables(self.variables_discharge_sim)
        self.var_time_obs, self.var_discharge_obs, self.var_wlevel_obs = self.define_file_variables(self.variables_obs)

        self.freq_discharge = 'H'
        self.periods_discharge_from = 72
        self.periods_discharge_to = 24
        self.file_time_discharge = self.define_file_time()

        self.folder_name_anc_sim = anc_dict[self.flag_discharge_data_sim][self.folder_name_tag]
        self.file_name_anc_sim = anc_dict[self.flag_discharge_data_sim][self.file_name_tag]
        self.folder_name_anc_obs = anc_dict[self.flag_discharge_data_obs][self.folder_name_tag]
        self.file_name_anc_obs = anc_dict[self.flag_discharge_data_obs][self.file_name_tag]

        self.file_path_anc_sim = self.define_file_ancillary(
            self.time_now, self.folder_name_anc_sim, self.file_name_anc_sim)

        self.file_path_anc_obs = self.define_file_ancillary(
            self.time_now, self.folder_name_anc_obs, self.file_name_anc_obs)

        self.flag_cleaning_anc_discharge_sim = flag_cleaning_anc_discharge_sim
        self.flag_cleaning_anc_discharge_obs = flag_cleaning_anc_discharge_obs
        # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to define file variable(s)
    def define_file_variables(self, variables_obj):
        var_time, var_discharge, var_water_level = None, None, None
        for var_key, var_name in variables_obj.items():
            if var_key == self.var_name_time:
                var_time = var_name
            elif var_key == self.var_name_discharge:
                var_discharge = var_name
            elif var_key == self.var_name_water_level:
                var_water_level = var_name
        return var_time, var_discharge, var_water_level
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to define time period
    def define_file_time(self):

        time_run = self.time_run

        time_day_start = time_run.replace(hour=0)
        time_day_end = time_run.replace(hour=23)

        time_period_from = pd.date_range(
            end=time_day_start, periods=self.periods_discharge_from, freq=self.freq_discharge)
        time_period_day = pd.date_range(
            start=time_day_start, end=time_day_end, freq=self.freq_discharge)
        time_period_to = pd.date_range(
            start=time_day_end, periods=self.periods_discharge_to, freq=self.freq_discharge)

        time_period = time_period_from.union(time_period_day).union(time_period_to)

        return time_period

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to define ancillary filename
    def define_file_ancillary(self, time, folder_name_raw, file_name_raw):

        alg_template_tags = self.alg_template_tags

        file_path_dict = {}
        for domain_name in self.domain_name_list:

            alg_template_values = {'domain_name': domain_name,
                                   'ancillary_sub_path_time_discharge': time,
                                   'ancillary_datetime_discharge': time}

            folder_name_def = fill_tags2string(folder_name_raw, alg_template_tags, alg_template_values)
            file_name_def = fill_tags2string(file_name_raw, alg_template_tags, alg_template_values)

            file_path_def = os.path.join(folder_name_def, file_name_def)

            file_path_dict[domain_name] = file_path_def

        return file_path_dict

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to define discharge filename
    def define_file_discharge(self, time, folder_name_raw, file_name_raw,
                              file_name_prefix=None, file_name_suffix=None, file_name_sep='_', file_name_elem=3,
                              file_sort_descending=True, extra_args=None):

        alg_template_tags = self.alg_template_tags

        section_name_obj = None
        time_period, time_rounding, time_frequency = None, None, None
        if extra_args is not None:
            if 'section_name_obj' in list(extra_args.keys()):
                section_name_obj = extra_args['section_name_obj']
            if 'time_period' in list(extra_args.keys()):
                time_period = extra_args['time_period']
            if 'time_rounding' in list(extra_args.keys()):
                time_rounding = extra_args['time_rounding']
            if 'time_frequency' in list(extra_args.keys()):
                time_frequency = extra_args['time_frequency']
            if 'section_description_obj' in list(extra_args.keys()):
                section_description_obj = extra_args['section_description_obj']

        if (time_rounding is not None) and (time_period is not None):
            time_range = pd.date_range(end=time, periods=time_period, freq=time_frequency)
            time_start = time_range[0].floor(time_rounding)
            time_end = time_range[-1]
        else:
            time_start = time
            time_end = time

        file_path_dict = {}
        for domain_name in self.domain_name_list:

            file_path_dict[domain_name] = {}

            section_name_list = None
            if section_name_obj is not None:
                if domain_name in list(section_name_obj.keys()):
                    section_name_list = section_name_obj[domain_name]

            alg_template_values = {'domain_name': domain_name,
                                   #"section_description": section_description_obj[domain_name],
                                   'source_sub_path_time_discharge_sim': time,
                                   'source_datetime_discharge_sim': time,
                                   'source_datetime_from_discharge_sim': '*',
                                   'source_datetime_to_discharge_sim': time_end,
                                   'source_sub_path_time_discharge_obs': time,
                                   'source_datetime_discharge_obs': time,
                                   'source_datetime_from_discharge_obs': '*',
                                   'source_datetime_to_discharge_obs': time_end,
                                   'ancillary_sub_path_time_discharge': time,
                                   'ancillary_datetime_discharge': time,
                                   'mask_name': '*',
                                   'scenario_discharge': '*'}

            if section_name_list is None:
                folder_name_def = fill_tags2string(folder_name_raw, alg_template_tags, alg_template_values)
                file_name_def = fill_tags2string(file_name_raw, alg_template_tags, alg_template_values)

                file_path_def = os.path.join(folder_name_def, file_name_def)

                section_path_found = glob.glob(file_path_def)
                if file_name_elem is not None:
                    section_path_obj = []
                    for section_path_step in section_path_found:
                        folder_name_step, file_name_step = os.path.split(section_path_step)
                        if (file_name_prefix is not None) and (file_name_suffix is None):
                            if file_name_step.startswith(file_name_prefix):

                                prefix_check = False
                                file_name_parts = file_name_step.split(file_name_sep)
                                file_prefix_parts = file_name_prefix.split(file_name_sep)
                                if file_name_parts.__len__() == file_name_elem:
                                    for prefix_id, prefix_step in enumerate(file_prefix_parts):
                                        if prefix_step == file_name_parts[prefix_id]:
                                            prefix_check = True
                                        else:
                                            prefix_check = False
                                            break
                                    if prefix_check:
                                        section_path_obj.append(section_path_step)

                        elif (file_name_suffix is not None) and (file_name_prefix is None):
                            if file_name_step.endswith(file_name_suffix):
                                section_path_obj.append(section_path_step)
                        else:
                            log_stream.error(' ===> Filter using "prefix" and "suffix" is not supported ')
                            raise NotImplementedError('Case not implemented yet')
                else:
                    section_path_obj = deepcopy(section_path_found)

                if not section_path_obj:

                    log_stream.error(' ===> Discharge simulated file are not available using the following ' +
                                     file_path_def + '. Try to use unfilled template string')

                    file_name_root = deepcopy(file_name_raw)
                    for template_key, template_value in alg_template_tags.items():
                        string_key = '{' + template_key + '}'

                        file_name_root = file_name_root.replace(string_key, '*')

                    file_part_start = file_name_root.split('*')[0]
                    file_part_end = file_name_root.split('*')[-1]

                    file_part_merge = ''
                    if file_part_start == file_part_end:
                        file_part_merge = file_part_start + '*'
                    elif file_part_start != file_part_end:
                        file_part_merge = file_part_start + '*' + file_part_end

                    log_stream.warning(' ===> Discharge simulated file are not available using the following ' +
                                       file_name_def + '. Try to use unfilled template string ' + file_part_merge)

                    file_path_def = os.path.join(folder_name_def, file_part_merge)
                    section_path_obj = glob.glob(file_path_def)

                section_path_obj.sort(reverse=file_sort_descending)
            else:
                section_path_obj = {}
                for section_name_step in section_name_list:

                    section_name_all, section_name_part1, section_name_part2 = None, None, None
                    if '_' in section_name_step:
                        section_name_parts = section_name_step.split('_')
                        if section_name_parts.__len__() == 2:
                            section_name_part1, section_name_part2 = section_name_parts[0], section_name_parts[1]
                        else:
                            log_stream.error(' ===> Section name must be defined by two elements')
                            raise RuntimeError('Define string using two elements with "_" separator')
                    else:
                        section_name_all = deepcopy(section_name_step)

                    alg_template_extra = {'section_name': section_name_all,
                                          'section_name_part1': section_name_part1,
                                          'section_name_part2': section_name_part2
                                          }

                    alg_template_values = {**alg_template_values, **alg_template_extra}

                    folder_name_def = fill_tags2string(folder_name_raw, alg_template_tags, alg_template_values)
                    file_name_def = fill_tags2string(file_name_raw, alg_template_tags, alg_template_values)

                    file_path_def = os.path.join(folder_name_def, file_name_def)

                    file_path_list = glob.glob(file_path_def)
                    file_path_list.sort(reverse=file_sort_descending)

                    section_path_obj[section_name_step] = file_path_list

            file_path_dict[domain_name] = section_path_obj

        return file_path_dict

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to wrap method(s)
    def organize_discharge(self):

        if self.scenario_type == 'simulated':
            section_collections = self.organize_discharge_sim()
        elif self.scenario_type == 'observed':
            section_collections = self.organize_discharge_obs()
        elif self.scenario_type == 'mixed':

            section_collections_sim = self.organize_discharge_sim()
            section_collections_obs = self.organize_discharge_obs()

            section_collections = self.organize_discharge_mixed(section_collections_obs, section_collections_sim)

        else:
            log_stream.error(' ===> Scenario type "' + self.scenario_type + '" is not expected')
            raise RuntimeError('Scenario type permitted flags are: [observed, simulated]')

        section_collections = self.filter_discharge(section_collections)

        return section_collections

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to filter the discharge collections
    def filter_discharge(self, domain_data_collections):

        time_run = self.time_run
        geo_data_collection = self.geo_data_collection

        log_stream.info(' ---> Filter discharge datasets [' + time_run.strftime(time_format_algorithm) + '] ... ')

        file_time_discharge = self.file_time_discharge
        scenario_boundary = self.scenario_boundary

        section_collection_filter = {}
        for domain_name_step in self.domain_name_list:

            log_stream.info(' ----> Domain "' + domain_name_step + '" ... ')

            section_data_collections = domain_data_collections[domain_name_step]
            section_geo_collection = geo_data_collection[domain_name_step]['section_data']

            log_stream.info(' -----> Collect datasets ... ')

            section_workspace_filter, section_workspace_valid = {}, {}
            time_first_list, time_last_list, idx_first_list, idx_last_list = [], [], [], []
            for section_key, section_data in section_geo_collection.items():

                section_description = section_data['section_description']

                log_stream.info(' ------> Section "' + section_description + '" ... ')

                section_workspace_valid[section_key] = {}
                if section_description in list(section_data_collections.keys()):
                    time_series_collections = section_data_collections[section_description]

                    if time_series_collections is not None:

                        time_first_valid = time_series_collections[self.var_name_discharge].first_valid_index()
                        idx_first_valid = time_series_collections[self.var_name_discharge].index.get_loc(time_first_valid)

                        time_last_valid = time_series_collections[self.var_name_discharge].last_valid_index()
                        idx_last_valid = time_series_collections[self.var_name_discharge].index.get_loc(time_last_valid)

                        time_first_list.append(time_first_valid)
                        time_last_list.append(time_last_valid)
                        idx_first_list.append(idx_first_valid)
                        idx_last_list.append(idx_last_valid)

                        time_series_attrs = {'time_first_valid': time_first_valid, 'time_last_valid': time_last_valid,
                                             'idx_first_valid': idx_first_valid, 'idx_last_valid': idx_last_valid}
                        time_series_collections_attrs = {**time_series_attrs, **time_series_collections.attrs}
                        time_series_collections.attrs = deepcopy(time_series_collections_attrs)

                        section_workspace_valid[section_description] = deepcopy(time_series_collections)

                        log_stream.info(' ------> Section "' + section_description +
                                        '" ... DONE')

                    else:
                        section_workspace_valid[section_description] = None

                        log_stream.info(' ------> Section "' + section_description +
                                        '" ... SKIPPED. Datasets are defined by NoneType ')

                else:
                    section_workspace_valid[section_description] = None

                    log_stream.info(' ------> Section "' + section_description +
                                    '" ... SKIPPED. Series are defined by NoneType ')

            log_stream.info(' -----> Collect datasets ... DONE')

            log_stream.info(' -----> Validate datasets ... ')
            if time_first_list and time_last_list:

                idx_first_select = max(idx_first_list)
                idx_last_select = min(idx_last_list)

                time_first_select = max(time_first_list)
                time_last_select = min(time_last_list)

                time_series_filter = pd.DataFrame(index=file_time_discharge)
                for section_key, section_data in section_geo_collection.items():

                    section_description = section_data['section_description']

                    log_stream.info(' -----> Section "' + section_description + '" ... ')

                    section_workspace_valid[section_key] = {}
                    if section_description in list(section_data_collections.keys()):
                        time_series_collections = section_data_collections[section_description]
                        if scenario_boundary is not None:
                            if scenario_boundary == 'both':
                                if time_series_collections is not None:
                                    time_series_selected = time_series_collections[time_first_select:time_last_select]
                                else:
                                    time_series_selected = None

                            elif scenario_boundary == 'right':
                                time_series_selected = time_series_collections[:time_last_select]
                            elif scenario_boundary == 'left':
                                time_series_selected = time_series_collections[time_first_select:]
                            else:
                                log_stream.error(' ===> Scenario boundary method "' + scenario_boundary +
                                                 '" is not supported')
                                raise NotImplementedError('Case not implemented yet')
                        else:
                            time_series_selected = deepcopy(time_series_collections)

                        if time_series_selected is not None:
                            time_series_filter[self.var_name_discharge] = time_series_selected[self.var_name_discharge]
                            time_series_filter[self.var_name_type] = time_series_selected[self.var_name_type]

                            time_series_filter_attrs = {**time_series_collections.attrs, **{'filter_stream': scenario_boundary}}
                            time_series_filter.attrs = deepcopy(time_series_filter_attrs)

                            section_workspace_filter[section_description] = deepcopy(time_series_filter)
                            log_stream.info(' -----> Section "' + section_description + '" ... DONE')
                        else:
                            section_workspace_filter[section_description] = None
                            log_stream.info(' ------> Section "' + section_description +
                                            '" ... SKIPPED. Series are defined by NoneType ')


                    else:
                        section_workspace_filter[section_description] = None

                        log_stream.info(' ------> Section "' + section_description +
                                        '" ... SKIPPED. Series are defined by NoneType ')

                section_collection_filter[domain_name_step] = section_workspace_filter
                log_stream.info(' -----> Validate datasets ... DONE')
            else:

                section_collection_filter[domain_name_step] = None
                log_stream.info(' -----> Validate datasets ... SKIPPED. All datasets are NoneType')

            log_stream.info(' ----> Domain "' + domain_name_step + '" ... DONE')

            log_stream.info(' ---> Filter discharge datasets [' + time_run.strftime(time_format_algorithm) +
                            '] ... DONE')

        return section_collection_filter
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to mix discharge collection(s)
    def organize_discharge_mixed(self, section_collection_obs, section_collection_sim):

        time_run = self.time_run
        geo_data_collection = self.geo_data_collection

        log_stream.info(' ---> Organize mixed discharge datasets [' + time_run.strftime(time_format_algorithm) +
                        '] ... ')

        section_collection_mix = {}
        for domain_name_step in self.domain_name_list:

            log_stream.info(' ----> Domain "' + domain_name_step + '" ... ')

            section_obs = section_collection_obs[domain_name_step]
            section_sim = section_collection_sim[domain_name_step]

            section_geo_collection = geo_data_collection[domain_name_step]['section_data']

            section_workspace_mix = {}
            for section_key, section_data in section_geo_collection.items():

                section_description = section_data['section_description']

                log_stream.info(' -----> Section "' + section_description + '" ... ')

                if section_description in list(section_obs.keys()):
                    time_series_obs = section_obs[section_description]
                    if time_series_obs is not None:
                        attrs_obs = time_series_obs.attrs
                        attr_links_obs = attrs_obs['link_stream']
                    else:
                        log_stream.warning(' ===> Section time-series is undefined in the observed datasets')
                        time_series_obs, attr_links_obs = None, None
                else:
                    log_stream.warning(' ===> Section time-series is not available in the observed datasets')
                    time_series_obs, attr_links_obs = None, None

                if section_description in list(section_sim.keys()):
                    time_series_sim = section_sim[section_description]
                    if time_series_sim is not None:
                        attrs_sim = time_series_sim.attrs
                        attr_links_sim = attrs_sim['link_stream']
                    else:
                        log_stream.warning(' ===> Section time-series is undefined in the simulated datasets')
                        time_series_sim, attr_links_sim = None, None
                else:
                    log_stream.warning(' ===> Section time-series is not available in the simulated datasets')
                    time_series_sim, attr_links_sim = None, None

                log_stream.info(' ------> Organize mixed datasets ... ')
                if attr_links_obs:

                    if (time_series_obs is not None) and (time_series_sim is not None):

                        log_stream.info(' -------> Case mixed with simulated values only ... ')
                        log_stream.info(' --------> Observed dset links == True')
                        log_stream.info(' --------> Observed dset data == defined :: Simulated dset data == defined')

                        time_period_mix = list(time_series_sim.index)
                        value_period_mix = list(time_series_sim[self.var_name_discharge].values)
                        type_period_mix = list(time_series_sim[self.var_name_type].values)

                        section_dframe_mix = create_obj_hydro(
                            time_period_mix, [self.var_name_discharge, self.var_name_type],
                            [value_period_mix, type_period_mix])

                        attrs_mix = {'link_stream': attr_links_sim, 'type_stream': 'simulated'}
                        section_dframe_mix.attrs = attrs_mix

                        log_stream.info(' -------> Case mixed with simulated values only ... DONE')

                    elif (time_series_obs is not None) and (time_series_sim is None):

                        log_stream.info(' -------> Case mixed with observed values only ... ')
                        log_stream.info(' --------> Observed dset links == True')
                        log_stream.info(' --------> Observed dset data == defined :: Simulated dset data == undefined')

                        time_period_obs = list(time_series_obs.index)
                        value_period_obs = list(time_series_obs[self.var_name_discharge].values)
                        type_period_obs = list(time_series_obs[self.var_name_type].values)

                        section_dframe_mix = create_obj_hydro(
                            time_period_obs, [self.var_name_discharge, self.var_name_type],
                            [value_period_obs, type_period_obs])

                        attrs_mix = {'link_stream': attr_links_obs, 'type_stream': 'observed'}
                        section_dframe_mix.attrs = deepcopy(attrs_mix)

                        log_stream.info(' -------> Case mixed with observed values only ... DONE')

                    else:

                        log_stream.info(' -------> Case mixed with undefined values ... ')
                        log_stream.info(' --------> Observed dset links == True')
                        log_stream.info(' --------> Observed dset data == undefined :: Simulated dset data == undefined')
                        log_stream.warning(' ===> All datasets are undefined or NoneType')
                        section_dframe_mix = None
                        log_stream.info(' -------> Case mixed with undefined values ... DONE')

                else:

                    if (time_series_obs is not None) and (time_series_sim is not None):

                        log_stream.info(' -------> Case mixed with observed/simulated values ... ')
                        log_stream.info(' --------> Observed dset links == False')
                        log_stream.info(' --------> Observed dset data == defined :: Simulated dset data == defined')

                        time_period_obs = list(time_series_obs.index)
                        value_period_obs = list(time_series_obs[self.var_name_discharge].values)
                        type_period_obs = list(time_series_obs[self.var_name_type].values)
                        time_period_sim = list(time_series_sim.index)
                        value_period_sim = list(time_series_sim[self.var_name_discharge].values)
                        type_period_sim = list(time_series_sim[self.var_name_type].values)

                        time_period_mix, value_period_mix, type_period_mix = [], [], []
                        for time_step_obs, time_step_sim, value_step_obs, value_step_sim, \
                            type_step_obs, type_step_sim in zip(time_period_obs, time_period_sim,
                                                                value_period_obs, value_period_sim,
                                                                type_period_obs, type_period_sim):

                            assert time_step_sim == time_step_obs

                            time_step_mix = time_step_sim

                            if np.isnan(value_step_obs) or (value_step_obs < 0.0):
                                if np.isnan(value_step_sim):
                                    value_step_mix = np.nan
                                    type_step_mix = np.nan
                                else:
                                    value_step_mix = value_step_sim
                                    type_step_mix = type_step_sim
                            else:
                                value_step_mix = value_step_obs
                                type_step_mix = type_step_obs

                            time_period_mix.append(time_step_mix)
                            value_period_mix.append(value_step_mix)
                            type_period_mix.append(type_step_mix)

                        section_dframe_mix = create_obj_hydro(
                            time_period_mix, [self.var_name_discharge, self.var_name_type],
                            [value_period_mix, type_period_mix])

                        attrs_mix = {'link_stream': [attr_links_obs, attr_links_sim], 'type_stream': 'observed,simulated'}
                        section_dframe_mix.attrs = deepcopy(attrs_mix)

                        log_stream.info(' -------> Case mixed with observed/simulated values ... DONE')

                    elif (time_series_obs is None) and (time_series_sim is not None):

                        log_stream.info(' -------> Case mixed with simulated values only ... ')
                        log_stream.info(' --------> Observed dset links == False')
                        log_stream.info(' --------> Observed dset data == undefined :: Simulated dset data == defined')

                        time_period_mix = list(time_series_sim.index)
                        value_period_mix = list(time_series_sim[self.var_name_discharge].values)
                        type_period_mix = list(time_series_sim[self.var_name_type].values)

                        section_dframe_mix = create_obj_hydro(
                            time_period_mix, [self.var_name_discharge, self.var_name_type],
                            [value_period_mix, type_period_mix])

                        attrs_mix = {'link_stream': attr_links_sim, 'type_stream': 'simulated'}
                        section_dframe_mix.attrs = attrs_mix

                        log_stream.info(' -------> Case mixed with simulated values only ... DONE')

                    else:

                        log_stream.info(' -------> Case mixed with undefined values ... ')
                        log_stream.info(' --------> Observed dset links == False')
                        log_stream.info(' --------> Observed dset data == undefined :: Simulated dset data == undefined')
                        log_stream.warning(' ===> All datasets are undefined or NoneType')
                        section_dframe_mix = None
                        log_stream.info(' -------> Case mixed with undefined values ... DONE')

                log_stream.info(' ------> Organize mixed datasets ... DONE')

                section_workspace_mix[section_description] = deepcopy(section_dframe_mix)

                log_stream.info(' -----> Section "' + section_description + '" ... DONE')

            section_collection_mix[domain_name_step] = section_workspace_mix

            log_stream.info(' ---> Organize mixed discharge datasets [' + time_run.strftime(time_format_algorithm) +
                            '] ... DONE')

        return section_collection_mix

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to organize simulated discharge
    def organize_discharge_sim(self, data_type='simulated'):

        time_run = self.time_run
        geo_data_collection = self.geo_data_collection

        log_stream.info(' ---> Organize simulated discharge datasets [' + time_run.strftime(time_format_algorithm) +
                        '] ... ')

        file_path_discharge_sim = self.file_path_discharge_sim
        file_path_ancillary_sim = self.file_path_anc_sim
        file_time_discharge = self.file_time_discharge

        section_data_collection = {}
        for domain_name_step in self.domain_name_list:

            log_stream.info(' ----> Domain "' + domain_name_step + '" ... ')

            file_path_discharge = file_path_discharge_sim[domain_name_step]
            file_path_ancillary = file_path_ancillary_sim[domain_name_step]

            if self.flag_cleaning_anc_discharge_sim:
                if os.path.exists(file_path_ancillary):
                    os.remove(file_path_ancillary)

            if not os.path.exists(file_path_ancillary):

                log_stream.info(' -----> Get datasets ... ')
                section_geo_collection = geo_data_collection[domain_name_step]['section_data']

                section_workspace = {}
                for section_key, section_data in section_geo_collection.items():

                    section_description = section_data['section_description']
                    section_name = section_data['name_point_outlet']

                    log_stream.info(' ------> Section "' + section_description + '" ... ')

                    if section_name in list(file_path_discharge.keys()):

                        file_path_data = file_path_discharge[section_name]
                        log_stream.info(' ------> File path(s): "' + ','.join(file_path_data) + '"')
                        driver_type = DriverType(
                            section_name, file_path_data,
                            file_time=file_time_discharge,
                            variables_names=self.variables_discharge_sim,
                            file_type=self.type_sim,
                            data_type=data_type,
                            method_data_occurrence=self.method_data_occurrence_sim,
                            method_data_filling=self.method_data_filling_sim)

                        section_dframe = driver_type.get_data()
                        if section_dframe is not None:
                            log_stream.info(' ------> Section "' + section_description + '" ... DONE')
                        else:
                            log_stream.info(' ------> Section "' + section_description +
                                            '" ... SKIPPED. Files not available')
                    else:
                        log_stream.info(' ------> Section "' + section_description +
                                        '" ... SKIPPED. Section not available')
                        section_dframe = None

                    section_workspace[section_description] = section_dframe

                log_stream.info(' -----> Get datasets ... DONE')

                log_stream.info(' -----> Analyze datasets ... ')
                section_workspace = analyze_obj_hydro(
                    section_workspace, file_source=data_type, file_method=self.method_data_analysis_sim)
                log_stream.info(' -----> Analyze datasets ... DONE')

                log_stream.info(' -----> Check datasets ... ')
                for section_key_step, section_dframe_step in section_workspace.items():

                    log_stream.info(' ------> Section "' + section_key_step + '" ... ')

                    if section_dframe_step is None:

                        log_stream.info(' -------> Get empty datasets ... ')

                        if self.method_data_null_sim == 'links':

                            section_fields_step = None
                            for section_tag_tmp, section_fields_tmp in section_geo_collection.items():
                                if section_fields_tmp['section_description'] == section_key_step:
                                    section_fields_step = section_fields_tmp.copy()
                                    break

                            if section_fields_step is not None:
                                section_alinks_up = section_fields_step['area_links']['upstream']
                                section_alinks_down = section_fields_step['area_links']['downstream']

                                if (section_alinks_up is not None) or (section_alinks_down is not None):

                                    section_alinks = None
                                    if section_alinks_up is not None:
                                        section_alinks = section_alinks_up
                                    elif section_alinks_down is not None:
                                        section_alinks = section_alinks_down

                                    section_values_pnt = None
                                    section_area_ratio_ref = None
                                    if section_alinks is not None:
                                        for section_key_alinks, section_values_alinks in section_alinks.items():
                                            if section_key_alinks in list(section_workspace.keys()):
                                                section_dframe_alinks = section_workspace[section_key_alinks]
                                                if section_dframe_alinks is not None:
                                                    section_idx = section_dframe_alinks.index
                                                    section_values_tmp = section_dframe_alinks[self.var_name_discharge].values
                                                    section_type_tmp = section_dframe_alinks[self.var_name_type].values
                                                else:
                                                    section_idx = None
                                                    section_values_tmp = None
                                                    section_type_tmp = None
                                            else:
                                                section_dframe_alinks = None
                                                section_idx = None
                                                section_values_tmp = None
                                                section_type_tmp = None

                                            if section_dframe_alinks is not None:
                                                section_area_ratio_pnt = section_values_alinks['area_ratio_pnt']
                                                section_area_ratio_ref = section_values_alinks['area_ratio_ref']

                                                section_values_tmp = section_values_tmp * section_area_ratio_pnt
                                                if section_values_pnt is None:
                                                    section_values_pnt = section_values_tmp.copy()
                                                else:
                                                    section_values_pnt = section_values_pnt + section_values_tmp
                                            else:
                                                section_values_pnt = None

                                        if (section_area_ratio_ref is not None) and (section_values_pnt is not None):
                                            section_values_ref = section_values_pnt * section_area_ratio_ref
                                            section_type_ref = section_type_tmp.copy()

                                            section_dframe_ref = create_obj_hydro(
                                                section_idx, [self.var_name_discharge, self.var_name_type],
                                                [section_values_ref, section_type_ref])
                                            section_dframe_ref.attrs = {'link_stream': True, 'type_stream': data_type}
                                        else:
                                            section_dframe_ref = None
                                    else:
                                        section_dframe_ref = None

                                elif (section_alinks_up is not None) and (section_alinks_down is not None):
                                    log_stream.error(
                                        ' ===> Section links for upstream and downstream conditions is not supported.')
                                    raise NotImplementedError('Case not implemented yet')
                                else:
                                    log_stream.error(
                                        ' ===> Section links for upstream and downstream conditions is not allowed.')
                                    raise RuntimeError('Check your upstream and downstream conditions')

                            else:
                                section_dframe_ref = None
                                log_stream.info(' -------> Get empty datasets ... FAILED. '
                                                'Datasets of the other points are empty.')

                            log_stream.info(
                                ' -------> Get empty datasets ... filled by using upstream and downstream conditions ')

                        else:
                            section_dframe_ref = None
                            log_stream.info(' -------> Get empty datasets ... SKIPPED. '
                                            'None of the methods to filling datasets was selected.')

                        section_workspace[section_key_step] = section_dframe_ref

                    else:
                         section_workspace[section_key_step].attrs = {'link_stream': False, 'type_stream': data_type}

                    log_stream.info(' ------> Section "' + section_key_step + '" ... DONE')

                log_stream.info(' -----> Check datasets ... DONE')

                flag_save_obj = True
                for section_key, section_data in section_workspace.items():
                    if section_data is None:
                        flag_save_obj = False
                        break

                if flag_save_obj:

                    folder_name_ancillary, file_name_ancillary = os.path.split(file_path_ancillary)
                    make_folder(folder_name_ancillary)

                    write_obj(file_path_ancillary, section_workspace)
                    log_stream.info(' ----> Domain "' + domain_name_step + '" ... DONE')
                else:
                    log_stream.info(' ----> Domain "' + domain_name_step +
                                    '" ... SKIPPED. All or some datasets are empty')

            else:

                section_workspace = read_obj(file_path_ancillary)

                log_stream.info(' ----> Domain "' + domain_name_step + '" ... SKIPPED. Data previously computed')

            section_data_collection[domain_name_step] = section_workspace

            log_stream.info(' ---> Organize simulated discharge datasets [' + time_run.strftime(time_format_algorithm) +
                            '] ... DONE')

        return section_data_collection

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to organize observed discharge
    def organize_discharge_obs(self, data_type='observed'):

        time_run = self.time_run
        geo_data_collection = self.geo_data_collection

        log_stream.info(' ---> Organize observed discharge datasets [' + time_run.strftime(time_format_algorithm) +
                        '] ... ')

        file_path_discharge_obs = self.file_path_discharge_obs
        file_path_ancillary_obs = self.file_path_anc_obs
        file_time_discharge = self.file_time_discharge

        section_data_collection = {}
        for domain_name_step in self.domain_name_list:

            log_stream.info(' ----> Domain "' + domain_name_step + '" ... ')

            file_path_discharge = file_path_discharge_obs[domain_name_step]
            file_path_ancillary = file_path_ancillary_obs[domain_name_step]

            if self.flag_cleaning_anc_discharge_obs:
                if os.path.exists(file_path_ancillary):
                    os.remove(file_path_ancillary)

            if not os.path.exists(file_path_ancillary):

                log_stream.info(' -----> Get datasets ... ')
                section_geo_collection = geo_data_collection[domain_name_step]['section_data']

                section_workspace = {}
                for section_key, section_data in section_geo_collection.items():

                    section_description = section_data['section_description']
                    section_name = section_data['name_point_obs']

                    log_stream.info(' ------> Section "' + section_description + '" ... ')

                    if section_name in list(file_path_discharge.keys()):

                        file_path_data = file_path_discharge[section_name]

                        driver_type = DriverType(
                            section_name, file_path_data,
                            file_time=file_time_discharge,
                            variables_names=self.variables_obs,
                            file_type=self.type_obs,
                            data_type=data_type,
                            method_data_occurrence=self.method_data_occurrence_obs,
                            method_data_filling=self.method_data_filling_obs)

                        section_dframe = driver_type.get_data()

                        log_stream.info(' ------> Section "' + section_description + '" ... DONE')
                    else:
                        log_stream.info(' ------> Section "' + section_description + '" ... SKIPPED. Datasets are empty')
                        section_dframe = None

                    section_workspace[section_description] = section_dframe

                log_stream.info(' -----> Get datasets ... DONE')

                log_stream.info(' -----> Analyze datasets ... ')
                section_workspace = analyze_obj_hydro(
                    section_workspace, file_source=data_type, file_method=self.method_data_analysis_obs)
                log_stream.info(' -----> Analyze datasets ... DONE')

                log_stream.info(' -----> Check datasets ... ')
                for section_key_step, section_dframe_step in section_workspace.items():

                    log_stream.info(' ------> Section "' + section_key_step + '" ... ')

                    if section_dframe_step is None:

                        log_stream.info(' -------> Get empty datasets ... ')

                        if self.method_data_null_obs == 'links':

                            section_fields_step = None
                            for section_tag_tmp, section_fields_tmp in section_geo_collection.items():
                                if section_fields_tmp['section_description'] == section_key_step:
                                    section_fields_step = section_fields_tmp.copy()
                                    break

                            if section_fields_step is not None:
                                section_alinks_up = section_fields_step['area_links']['upstream']
                                section_alinks_down = section_fields_step['area_links']['downstream']

                                if (section_alinks_up is not None) or (section_alinks_down is not None):

                                    section_alinks = None
                                    if section_alinks_up is not None:
                                        section_alinks = section_alinks_up
                                    elif section_alinks_down is not None:
                                        section_alinks = section_alinks_down

                                    section_values_pnt = None
                                    section_area_ratio_ref = None
                                    if section_alinks is not None:
                                        for section_key_alinks, section_values_alinks in section_alinks.items():
                                            if section_key_alinks in list(section_workspace.keys()):
                                                section_dframe_alinks = section_workspace[section_key_alinks]
                                                if section_dframe_alinks is not None:
                                                    section_idx = section_dframe_alinks.index
                                                    section_values_tmp = section_dframe_alinks[self.var_name_discharge].values
                                                    section_type_tmp = section_dframe_alinks[self.var_name_type].values
                                                else:
                                                    section_idx = None
                                                    section_values_tmp = None
                                                    section_type_tmp = None
                                            else:
                                                section_dframe_alinks = None
                                                section_idx = None
                                                section_values_tmp = None
                                                section_type_tmp = None

                                            if section_dframe_alinks is not None:
                                                section_area_ratio_pnt = section_values_alinks['area_ratio_pnt']
                                                section_area_ratio_ref = section_values_alinks['area_ratio_ref']

                                                section_values_tmp = section_values_tmp * section_area_ratio_pnt
                                                if section_values_pnt is None:
                                                    section_values_pnt = section_values_tmp.copy()
                                                else:
                                                    section_values_pnt = section_values_pnt + section_values_tmp
                                            else:
                                                section_values_pnt = None

                                        if (section_area_ratio_ref is not None) and (section_values_pnt is not None):
                                            section_values_ref = section_values_pnt * section_area_ratio_ref
                                            section_type_ref = section_type_tmp.copy()

                                            section_dframe_ref = create_obj_hydro(
                                                section_idx, [self.var_name_discharge, self.var_name_type],
                                                [section_values_ref, section_type_ref])
                                            section_dframe_ref.attrs = {'link_stream': True, 'type_stream': data_type}

                                            log_stream.info(
                                                ' -------> Get empty datasets ... DONE. '
                                                'Datasets filled by using upstream and downstream conditions')
                                        else:
                                            section_dframe_ref = None
                                    else:
                                        section_dframe_ref = None

                                elif (section_alinks_up is not None) and (section_alinks_down is not None):
                                    log_stream.error(
                                        ' ===> Section links for upstream and downstream conditions is not supported.')
                                    raise NotImplementedError('Case not implemented yet')
                                else:
                                    log_stream.error(
                                        ' ===> Section links for upstream and downstream conditions is not allowed.')
                                    raise RuntimeError('Check your upstream and downstream conditions')

                            else:
                                section_dframe_ref = None
                                log_stream.info(' -------> Get empty datasets ... DONE. '
                                                'Datasets are not filled by using a methods in the settings file')

                        else:
                            section_dframe_ref = None
                            log_stream.info(' -------> Get empty datasets ... SKIPPED. '
                                            'None of the methods to filling datasets was selected.')

                        section_workspace[section_key_step] = section_dframe_ref

                    else:
                        section_workspace[section_key_step].attrs = {'link_stream': False, 'type_stream': data_type}

                    log_stream.info(' ------> Section "' + section_key_step + '" ... DONE')

                log_stream.info(' -----> Check datasets ... DONE')

                # Save datasets
                flag_save_obj = True
                for section_key, section_data in section_workspace.items():
                    if section_data is None:
                        flag_save_obj = False
                        break

                if flag_save_obj:

                    folder_name_ancillary, file_name_ancillary = os.path.split(file_path_ancillary)
                    make_folder(folder_name_ancillary)

                    write_obj(file_path_ancillary, section_workspace)
                    log_stream.info(' ----> Domain "' + domain_name_step + '" ... DONE')
                else:
                    log_stream.info(' ----> Domain "' + domain_name_step +
                                    '" ... SKIPPED. All or some datasets are empty')

            else:

                section_workspace = read_obj(file_path_ancillary)
                log_stream.info(' ----> Domain "' + domain_name_step + '" ... SKIPPED. Data previously computed')

            section_data_collection[domain_name_step] = section_workspace

            log_stream.info(' ---> Organize observed discharge datasets [' + time_run.strftime(time_format_algorithm) +
                            '] ... DONE')

        return section_data_collection

    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
