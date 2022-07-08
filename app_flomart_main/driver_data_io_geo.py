"""
Class Features

Name:          driver_data_io_geo
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20200515'
Version:       '1.0.0'
"""

######################################################################################
# Library
import logging
import os
import numpy as np

from copy import deepcopy

from lib_utils_geo import read_file_geo, check_epsg_code
from lib_utils_hydraulic import read_file_hydraulic

from lib_utils_io import read_obj, write_obj
from lib_utils_system import fill_tags2string, make_folder
from lib_utils_generic import convert_array_2_list
from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
######################################################################################


# -------------------------------------------------------------------------------------
# Class DriverGeo
class DriverGeo:

    # -------------------------------------------------------------------------------------
    # Initialize class
    def __init__(self, src_dict, dst_dict,
                 alg_ancillary=None, alg_template_tags=None,
                 flag_geo_data='geo_data',
                 flag_telemac_data='telemac_data', flag_hazard_data='hazard_data',
                 flag_geo_data_in='geo_data', flag_geo_data_out='geo_data',
                 flag_cleaning_geo=True):

        self.flag_geo_data = flag_geo_data
        self.flag_telemac_data = flag_telemac_data
        self.flag_hazard_data = flag_hazard_data
        self.flag_geo_data_in = flag_geo_data_in
        self.flag_geo_data_out = flag_geo_data_out

        self.alg_ancillary = alg_ancillary

        self.alg_template_tags = alg_template_tags
        self.file_name_tag = 'file_name'
        self.folder_name_tag = 'folder_name'

        self.domain_name_list = self.alg_ancillary['domain_name']

        self.folder_name_geo_in_generic = src_dict[self.flag_geo_data_in]['generic'][self.folder_name_tag]
        self.file_name_geo_in_generic = src_dict[self.flag_geo_data_in]['generic'][self.file_name_tag]
        self.folder_name_geo_in_hydraulic = src_dict[self.flag_geo_data_in]['hydraulic'][self.folder_name_tag]
        self.file_name_geo_in_hydraulic = src_dict[self.flag_geo_data_in]['hydraulic'][self.file_name_tag]

        self.folder_name_geo_out = dst_dict[self.flag_geo_data_out][self.folder_name_tag]
        self.file_name_geo_out = dst_dict[self.flag_geo_data_out][self.file_name_tag]

        self.collections_geo = self.read_geo_info_generic()
        self.collections_hydraulic = self.read_geo_info_hydraulic()

        self.flag_cleaning_geo = flag_cleaning_geo

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to read geographical hydraulic information
    def read_geo_info_hydraulic(self):

        domain_name_list = self.domain_name_list
        template_tags = self.alg_template_tags

        folder_name_raw, file_name_raw = self.folder_name_geo_in_hydraulic, self.file_name_geo_in_hydraulic

        hydraulic_collections = {}
        for domain_name in domain_name_list:

            template_values = {'domain_name': domain_name}

            folder_name_def = fill_tags2string(folder_name_raw, template_tags, template_values)
            file_name_def = fill_tags2string(file_name_raw, template_tags, template_values)

            if os.path.exists(os.path.join(folder_name_def, file_name_def)):
                if file_name_def.endswith('json'):
                    hydraulic_data = read_file_hydraulic(os.path.join(folder_name_def, file_name_def))
                    hydraulic_collections[domain_name] = hydraulic_data
                else:
                    log_stream.error(' ===> Hydraulic section file ' + file_name_def + ' format is not supported')
                    raise NotImplementedError('Case not implemented yet')
            else:
                log_stream.error(' ===> Hydraulic section file ' + file_name_def + ' not available')
                raise IOError('Check your configuration file')

        return hydraulic_collections

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to read geographical generic information
    def read_geo_info_generic(self):

        domain_name_list = self.domain_name_list

        template_tags = self.alg_template_tags

        folder_name_raw, file_name_raw = self.folder_name_geo_in_generic, self.file_name_geo_in_generic

        geo_collections = {}
        for domain_name in domain_name_list:

            template_values = {'domain_name': domain_name}

            folder_name_def = fill_tags2string(folder_name_raw, template_tags, template_values)
            file_name_def = fill_tags2string(file_name_raw, template_tags, template_values)

            if os.path.exists(os.path.join(folder_name_def, file_name_def)):
                if file_name_def.endswith('mat'):
                    geo_data = read_file_geo(os.path.join(folder_name_def, file_name_def))
                    geo_collections[domain_name] = geo_data
                else:
                    log_stream.error(' ===> Geographical section file ' + file_name_def + ' format is not supported')
                    raise NotImplementedError('Case not implemented yet')
            else:
                log_stream.error(' ===> Geographical reference file ' + file_name_def + ' not available')
                raise IOError('Check your configuration file')
        return geo_collections

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to define domain collection
    def define_domain_collection(self, domain_name, folder_name, file_name):

        template_tags = self.alg_template_tags
        template_values = {'domain_name': domain_name}

        folder_name = fill_tags2string(folder_name, template_tags, template_values)
        file_name = fill_tags2string(file_name, template_tags, template_values)
        file_path = os.path.join(folder_name, file_name)

        if self.flag_cleaning_geo:
            if os.path.exists(file_path):
                os.remove(file_path)

        return file_path
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to define upstream and downstream link(s)
    @staticmethod
    def organize_section_links(section_info):

        section_drainage_area_tags = ['section_drainage_area', 'drainage_area']

        for section_tag_ref, section_fields_ref in section_info.items():
            section_name_upstream_ref = section_fields_ref['name_point_upstream']
            section_name_downstream_ref = section_fields_ref['name_point_downstream']
            section_name_ref = section_fields_ref['section_description']

            section_drainage_area_ref = None
            for section_drainage_area_tag in section_drainage_area_tags:
                if section_drainage_area_tag in list(section_fields_ref.keys()):
                    section_drainage_area_ref = section_fields_ref[section_drainage_area_tag]
                    break

            section_obj_summary = {'downstream': {}, 'upstream': {}}
            if section_name_downstream_ref is not None:

                section_drainage_area_tot = 0
                for section_name_step in section_name_downstream_ref:

                    section_obj_summary['downstream'][section_name_step] = {}
                    section_fields_step = None

                    for section_tag_tmp, section_fields_tmp in section_info.items():
                        if section_fields_tmp['section_description'] == section_name_step:
                            section_fields_step = section_fields_tmp.copy()
                            break

                    section_drainage_area_step = None
                    for section_drainage_area_tag in section_drainage_area_tags:
                        if section_fields_step is not None:
                            if section_drainage_area_tag in list(section_fields_step.keys()):
                                section_drainage_area_step = section_fields_step[section_drainage_area_tag]
                                break

                    if section_drainage_area_step is not None:
                        section_drainage_area_tot = section_drainage_area_tot + section_drainage_area_step
                        section_obj_summary['downstream'][section_name_step]['area_partial_pnt'] = section_drainage_area_step
                    else:
                        section_obj_summary['downstream'][section_name_step]['area_partial_pnt'] = None

                for section_name_step in section_name_downstream_ref:
                    section_drainage_area_step = section_obj_summary['downstream'][section_name_step]['area_partial_pnt']
                    if (section_drainage_area_step is not None) and (section_drainage_area_ref is not None):
                        section_ratio_pnt_step = section_drainage_area_step / section_drainage_area_tot
                        section_ratio_ref_step = section_drainage_area_ref / section_drainage_area_tot
                        section_obj_summary['downstream'][section_name_step]['area_ratio_pnt'] = section_ratio_pnt_step
                        section_obj_summary['downstream'][section_name_step]['area_ratio_ref'] = section_ratio_ref_step
                        section_obj_summary['downstream'][section_name_step]['area_total_ref'] = section_drainage_area_ref
                        section_obj_summary['downstream'][section_name_step]['area_total_pnt'] = section_drainage_area_tot
                    else:
                        section_obj_summary['downstream'][section_name_step]['area_ratio_pnt'] = None
                        section_obj_summary['downstream'][section_name_step]['area_ratio_ref'] = None
                        section_obj_summary['downstream'][section_name_step]['area_total_ref'] = None
                        section_obj_summary['downstream'][section_name_step]['area_total_pnt'] = None
            else:
                section_obj_summary['downstream'] = None

            if section_name_upstream_ref is not None:

                section_drainage_area_tot = 0
                for section_name_step in section_name_upstream_ref:

                    section_obj_summary['upstream'][section_name_step] = {}
                    section_fields_step = None

                    for section_tag_tmp, section_fields_tmp in section_info.items():
                        if section_fields_tmp['section_description'] == section_name_step:
                            section_fields_step = section_fields_tmp.copy()
                            break

                    section_drainage_area_step = None
                    for section_drainage_area_tag in section_drainage_area_tags:
                        if section_drainage_area_tag in list(section_fields_step.keys()):
                            section_drainage_area_step = section_fields_step[section_drainage_area_tag]
                            break

                    section_drainage_area_tot = section_drainage_area_tot + section_drainage_area_step
                    section_obj_summary['upstream'][section_name_step]['area_partial_pnt'] = section_drainage_area_step

                for section_name_step in section_name_upstream_ref:
                    section_drainage_area_step = section_obj_summary['upstream'][section_name_step]['area_partial_pnt']
                    if (section_drainage_area_step is not None) and (section_drainage_area_ref is not None):
                        section_ratio_pnt_step = section_drainage_area_step / section_drainage_area_tot
                        section_ratio_ref_step = section_drainage_area_ref / section_drainage_area_tot
                        section_obj_summary['upstream'][section_name_step]['area_ratio_pnt'] = section_ratio_pnt_step
                        section_obj_summary['upstream'][section_name_step]['area_ratio_ref'] = section_ratio_ref_step
                        section_obj_summary['upstream'][section_name_step]['area_total_ref'] = section_drainage_area_ref
                        section_obj_summary['upstream'][section_name_step]['area_total_pnt'] = section_drainage_area_tot
                    else:
                        section_obj_summary['upstream'][section_name_step]['area_ratio_pnt'] = None
                        section_obj_summary['upstream'][section_name_step]['area_ratio_ref'] = None
                        section_obj_summary['upstream'][section_name_step]['area_total_ref'] = None
                        section_obj_summary['upstream'][section_name_step]['area_total_pnt'] = None
            else:
                section_obj_summary['upstream'] = None

            section_info[section_tag_ref]['area_links'] = section_obj_summary

        return section_info

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to get section geo information
    @staticmethod
    def organize_section_geo(geo_data, geo_fields=None, section_key='section_{:}'):

        if geo_fields is None:
            geo_fields = ['section_description', 'section_id', 'section_name',
                          'section_drainage_area', 'section_altitude', 'section_discharge_idx']

        section_dset, section_collections, section_n = {}, {}, None
        for section_field in geo_fields:

            section_list = None
            if section_field in list(geo_data.keys()):
                section_value = geo_data[section_field]
                section_list = convert_array_2_list(section_value)
                section_dset[section_field] = section_list
            else:
                if section_field == 'section_altitude':
                    section_dset[section_field] = None
                elif section_field == 'section_drainage_area':
                    section_dset[section_field] = None
                else:
                    log_stream.error(' ===> Geo field "' + section_field + '" must be defined')
                    raise IOError('Field not found')

            if section_list is not None:
                section_n = section_list.__len__()

        if section_n is not None:

            for field_key, field_values in section_dset.items():

                for n in range(0, section_n):

                    if field_values is None:
                        field_values = [-9999] * section_n

                    field_value = field_values[n]

                    section_tag = section_key.format(n)
                    if section_tag not in list(section_collections.keys()):
                        section_collections[section_tag] = {}
                    section_collections[section_tag][field_key] = {}
                    section_collections[section_tag][field_key] = field_value

        return section_collections

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to get section hydraulic information
    @staticmethod
    def organize_section_hydraulic(section_data, hydraulic_data, hydraulic_key_expected=None):

        if hydraulic_key_expected is None:
            hydraulic_key_expected = ['name_point_outlet', 'name_point_downstream',
                                      'name_point_upstream', 'name_point_obs']

        for section_key, section_fields in section_data.items():
            section_description = section_fields['section_description']

            hydraulic_ws = None
            for hydraulic_id, hydraulic_fields in hydraulic_data.items():
                hydraulic_description = hydraulic_fields['description']
                if hydraulic_description == section_description:
                    hydraulic_ws = deepcopy(hydraulic_fields)
                    break

            if hydraulic_ws is not None:
                for hydraulic_key, hydraulic_value in hydraulic_ws.items():
                    if hydraulic_key not in list(section_fields.keys()):
                        if hydraulic_key in hydraulic_key_expected:
                            section_fields[hydraulic_key] = hydraulic_value
                        if (hydraulic_key == 'drainage_area') and (section_fields['section_drainage_area'] == -9999):
                            section_fields['section_drainage_area'] = hydraulic_value

                section_data[section_key] = section_fields

        return section_data

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to get map geo information
    @staticmethod
    def organize_map_geo(geo_data, geo_fields=None):

        if geo_fields is None:
            geo_fields = ['area_reference_id', 'area_mask', 'area_mask_extended',
                          'area_reference_geo_x', 'area_reference_geo_y',
                          'coord_bottom_y', 'coord_top_y', 'coord_left_x', 'coord_right_x',
                          'area_epsg_code']

        map_collections = {}
        for geo_field in geo_fields:

            if geo_field in list(geo_data.keys()):
                geo_value = geo_data[geo_field]

                if geo_value.ndim == 1:
                    if geo_value.shape[0] == 1:
                        geo_value = geo_value[0]
                        if isinstance(geo_value, str):
                            geo_value = str(geo_value)
                        else:
                            geo_value = float(geo_value)

                elif geo_value.ndim == 2:
                    if geo_value.shape[0] == 1 and geo_value.shape[1] == 1:
                        geo_value = float(geo_value[0][0])
                else:
                    log_stream.error(' ===> Geo field "' + geo_field + '" format is not supported')
                    raise NotImplementedError('Case not implemented yet')

                map_collections[geo_field] = geo_value
            else:

                if geo_field == 'area_reference_geo_x':
                    log_stream.warning(' ===> Geo field "' + geo_field + '" will be defined by default.')
                elif geo_field == 'area_reference_geo_y':
                    log_stream.warning(' ===> Geo field "' + geo_field + '" will be defined by default.')
                else:
                    log_stream.error(' ===> Geo field "' + geo_field + '" must be defined')
                    raise IOError('Field not found')

        if ('area_reference_geo_x' not in list(map_collections.items())) or \
                ('area_reference_geo_y' not in list(map_collections.items())):

            map_rows = map_collections['area_reference_id'].shape[0]
            map_cols = map_collections['area_reference_id'].shape[1]

            center_left_utm = map_collections['coord_left_x']
            center_right_utm = map_collections['coord_right_x']
            center_bottom_utm = map_collections['coord_bottom_y']
            center_top_utm = map_collections['coord_top_y']

            map_array_x = np.linspace(center_left_utm, center_right_utm, num=map_cols, endpoint=True)
            map_array_y = np.linspace(center_bottom_utm, center_top_utm, num=map_rows, endpoint=True)

            map_grid_x, map_grid_y = np.meshgrid(map_array_x, map_array_y)

            map_collections['area_reference_geo_x'] = map_grid_x
            map_collections['area_reference_geo_y'] = map_grid_y

        if 'area_epsg_code' in list(map_collections.keys()):
            epsg_code_unchecked = map_collections['area_epsg_code']
            epsg_code_checked = check_epsg_code(epsg_code_unchecked)
            map_collections['area_epsg_code'] = epsg_code_checked

        else:
            log_stream.error(' ===> Geo field "area_epsg_code" must be defined')
            raise IOError('Field not found')

        return map_collections

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to organize geographical data
    def organize_geo(self):

        log_stream.info(' ---> Organize geographical datasets ... ')

        geo_collection = {}
        for domain_name in self.domain_name_list:

            log_stream.info(' ----> Domain "' + domain_name + '" ... ')

            file_path_collections = self.define_domain_collection(
                domain_name, self.folder_name_geo_out, self.file_name_geo_out)
            geo_data = self.collections_geo[domain_name]
            hydraulic_data = self.collections_hydraulic[domain_name]

            if not os.path.exists(file_path_collections):

                section_data = self.organize_section_geo(geo_data)
                section_data = self.organize_section_hydraulic(section_data, hydraulic_data)
                section_data = self.organize_section_links(section_data)

                map_data = self.organize_map_geo(geo_data)

                obj_data = {'section_data': section_data, 'map_data': map_data}

                folder_name_collections, file_name_collections = os.path.split(file_path_collections)
                make_folder(folder_name_collections)
                write_obj(file_path_collections, obj_data)

                log_stream.info(' ----> Domain "' + domain_name + '" ... DONE')

            else:

                log_stream.info(' ----> Domain "' + domain_name + '" ... SKIPPED. Datasets previously computed.')
                obj_data = read_obj(file_path_collections)

            geo_collection[domain_name] = obj_data

        log_stream.info(' ---> Organize geographical datasets ... DONE')

        return geo_collection
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
