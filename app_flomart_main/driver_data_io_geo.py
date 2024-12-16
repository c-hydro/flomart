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

from lib_utils_section import join_name_and_alias
from lib_utils_geo import read_file_geo, check_epsg_code, write_geo_area
from lib_utils_hydraulic import read_file_hydraulic, map_description_section_2_hydraulic
from lib_utils_hydrology_ascii import read_file_hydro_section
from lib_utils_hydrology_generic import map_description_hydro_2_hydraulic

from lib_utils_tr import get_tr_params, get_tr_file
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
                 flag_hydro_data='hydro_data', flag_hazard_data='hazard_data',
                 flag_geo_data_in='geo_data',
                 flag_geo_data_out='geo_data', flag_geo_area_out='geo_area',
                 flag_geo_memory_area_out='geo_memory_area', flag_geo_memory_idx_out='geo_memory_idx',
                 flag_cleaning_geo=True):

        # define flags
        self.flag_geo_data = flag_geo_data
        self.flag_hydro_data = flag_hydro_data
        self.flag_hazard_data = flag_hazard_data

        self.flag_geo_data_in = flag_geo_data_in
        self.flag_geo_data_out = flag_geo_data_out
        self.flag_geo_area_out = flag_geo_area_out
        self.flag_geo_memory_area_out = flag_geo_memory_area_out
        self.flag_geo_memory_idx_out = flag_geo_memory_idx_out

        # set ancillary information
        self.alg_ancillary = alg_ancillary

        # set template tags and information tags
        self.alg_template_tags = alg_template_tags
        self.file_name_tag, self.folder_name_tag, self.file_type_tag = 'file_name', 'folder_name', 'file_type'

        # ancillary settings
        self.domain_name_list = self.alg_ancillary['domain_name']
        if 'tr_method' in list(self.alg_ancillary.keys()):
            self.method_tr = self.alg_ancillary['tr_method']
        else:
            self.method_tr = 'method_regional'

        if 'memory_saver' in list(self.alg_ancillary.keys()):
            self.memory_saver = self.alg_ancillary['memory_saver']
        else:
            self.memory_saver = False

        # define input datasets
        self.folder_name_geo_in_generic = src_dict[self.flag_geo_data_in]['generic'][self.folder_name_tag]
        self.file_name_geo_in_generic = src_dict[self.flag_geo_data_in]['generic'][self.file_name_tag]
        self.folder_name_geo_in_hydraulic = src_dict[self.flag_geo_data_in]['hydraulic'][self.folder_name_tag]
        self.file_name_geo_in_hydraulic = src_dict[self.flag_geo_data_in]['hydraulic'][self.file_name_tag]

        self.file_name_hydro_sections, self.folder_name_hydro_sections, self.file_type_hydro_sections = None, None, None
        self.file_name_hydro_qt, self.folder_name_hydro_qt, self.file_type_hydro_qt = None, None, None
        if self.flag_hydro_data in list(src_dict.keys()):
            self.file_name_hydro_sections = src_dict[self.flag_hydro_data]['sections'][self.file_name_tag]
            self.folder_name_hydro_sections = src_dict[self.flag_hydro_data]['sections'][self.folder_name_tag]
            self.file_type_hydro_sections = src_dict[self.flag_hydro_data]['sections'][self.file_type_tag]
            self.file_name_hydro_qt = src_dict[self.flag_hydro_data]['q_t'][self.file_name_tag]
            self.folder_name_hydro_qt = src_dict[self.flag_hydro_data]['q_t'][self.folder_name_tag]
            self.file_type_hydro_qt = src_dict[self.flag_hydro_data]['q_t'][self.file_type_tag]

        self.file_name_hazard = src_dict[self.flag_hazard_data][self.file_name_tag]
        self.folder_name_hazard = src_dict[self.flag_hazard_data][self.folder_name_tag]
        if 'file_type' in list(src_dict[self.flag_hazard_data].keys()):
            self.file_type_hazard = src_dict[self.flag_hazard_data][self.file_type_tag]
        else:
            self.file_type_hazard = 'tiff'

        # define outcome datasets
        self.folder_name_geo_data_out = dst_dict[self.flag_geo_data_out][self.folder_name_tag]
        self.file_name_geo_data_out = dst_dict[self.flag_geo_data_out][self.file_name_tag]

        if self.flag_geo_area_out in list(dst_dict.keys()):
            self.folder_name_geo_area_out = dst_dict[self.flag_geo_area_out][self.folder_name_tag]
            self.file_name_geo_area_out = dst_dict[self.flag_geo_area_out][self.file_name_tag]
        else:
            self.folder_name_geo_area_out = deepcopy(self.folder_name_geo_data_out)
            self.file_name_geo_area_out = 'area_{domain_name}.tiff'

        if self.flag_geo_memory_area_out in list(dst_dict.keys()):
            self.folder_name_geo_memory_area_out = dst_dict[self.flag_geo_memory_area_out][self.folder_name_tag]
            self.file_name_geo_memory_area_out = dst_dict[self.flag_geo_memory_area_out][self.file_name_tag]
        else:
            self.folder_name_geo_memory_area_out = deepcopy(self.folder_name_geo_data_out)
            self.file_name_geo_memory_area_out = 'memory_{domain_name}_data.workspace'

        if self.flag_geo_memory_idx_out in list(dst_dict.keys()):
            self.folder_name_geo_memory_idx_out = dst_dict[self.flag_geo_memory_idx_out][self.folder_name_tag]
            self.file_name_geo_memory_idx_out = dst_dict[self.flag_geo_memory_idx_out][self.file_name_tag]
        else:
            self.folder_name_geo_memory_idx_out = deepcopy(self.folder_name_geo_data_out)
            self.file_name_geo_memory_idx_out = 'memory_{domain_name}_idx_{section_id}.workspace'

        self.flag_cleaning_geo = flag_cleaning_geo

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to read hydrological information
    def read_hydro_info(self, domain_name):

        log_stream.info(' -----> Read hydrological information ... ')

        template_tags = self.alg_template_tags

        if (self.file_name_hydro_sections is not None) and (self.folder_name_hydro_sections is not None):

            folder_name_raw, file_name_raw = self.folder_name_hydro_sections, self.file_name_hydro_sections

            template_values = {'domain_name': domain_name}
            folder_name_def = fill_tags2string(folder_name_raw, template_tags, template_values)
            if isinstance(folder_name_def, tuple):
                folder_name_def = folder_name_def[0]
            file_name_def = fill_tags2string(file_name_raw, template_tags, template_values)
            if isinstance(file_name_def, tuple):
                file_name_def = file_name_def[0]

            if os.path.exists(os.path.join(folder_name_def, file_name_def)):
                if file_name_def.endswith('txt'):
                    hydrological_data = read_file_hydro_section(os.path.join(folder_name_def, file_name_def))
                else:
                    log_stream.error(' ===> Hydrological section file ' + file_name_def + ' format is not supported')
                    raise NotImplementedError('Case not implemented yet')
            else:
                log_stream.error(' ===> Hydrological section file ' + file_name_def + ' not available')
                raise IOError('Check your configuration file')

            log_stream.info(' -----> Read hydrological information ... DONE')

        else:
            log_stream.info(' -----> Read hydrological information ... DATA NOT DEFINED')
            hydrological_data = None

        return hydrological_data

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to read geographical hydraulic information
    def read_geo_info_hydraulic(self, domain_name):

        log_stream.info(' -----> Read geographical hydraulic information ... ')

        template_tags = self.alg_template_tags

        folder_name_raw, file_name_raw = self.folder_name_geo_in_hydraulic, self.file_name_geo_in_hydraulic

        template_values = {'domain_name': domain_name}
        folder_name_def = fill_tags2string(folder_name_raw, template_tags, template_values)
        if isinstance(folder_name_def, tuple):
            folder_name_def = folder_name_def[0]
        file_name_def = fill_tags2string(file_name_raw, template_tags, template_values)
        if isinstance(file_name_def, tuple):
            file_name_def = file_name_def[0]

        if os.path.exists(os.path.join(folder_name_def, file_name_def)):
            if file_name_def.endswith('json'):
                hydraulic_data = read_file_hydraulic(os.path.join(folder_name_def, file_name_def))
            else:
                log_stream.error(' ===> Hydraulic section file ' + file_name_def + ' format is not supported')
                raise NotImplementedError('Case not implemented yet')
        else:
            log_stream.error(' ===> Hydraulic section file ' + file_name_def + ' not available')
            raise IOError('Check your configuration file')

        log_stream.info(' -----> Read geographical hydraulic information ... DONE')

        return hydraulic_data

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to read geographical generic information
    def read_geo_info_generic(self, domain_name):

        log_stream.info(' -----> Read geographical generic information ... ')

        template_tags = self.alg_template_tags

        folder_name_raw, file_name_raw = self.folder_name_geo_in_generic, self.file_name_geo_in_generic

        template_values = {'domain_name': domain_name}
        folder_name_def = fill_tags2string(folder_name_raw, template_tags, template_values)
        if isinstance(folder_name_def, tuple):
            folder_name_def = folder_name_def[0]
        file_name_def = fill_tags2string(file_name_raw, template_tags, template_values)
        if isinstance(file_name_def, tuple):
            file_name_def = file_name_def[0]

        if os.path.exists(os.path.join(folder_name_def, file_name_def)):
            if file_name_def.endswith('mat'):
                geo_data = read_file_geo(os.path.join(folder_name_def, file_name_def))
            else:
                log_stream.error(' ===> Geographical section file ' + file_name_def + ' format is not supported')
                raise NotImplementedError('Case not implemented yet')
        else:
            log_stream.error(' ===> Geographical reference file ' + file_name_def + ' not available')
            raise IOError('Check your configuration file')

        log_stream.info(' -----> Read geographical generic information ... DONE')

        return geo_data

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to define domain collection
    def define_domain_collection(self, domain_name, folder_name, file_name, section_description=None):

        template_tags = self.alg_template_tags
        if section_description is not None:
            template_values = {'domain_name': domain_name, 'section_description': str(section_description)}
        else:
            template_values = {'domain_name': domain_name}

        folder_name = fill_tags2string(folder_name, template_tags, template_values)
        if isinstance(folder_name, tuple):
            folder_name = folder_name[0]
        file_name = fill_tags2string(file_name, template_tags, template_values)
        if isinstance(file_name, tuple):
            file_name = file_name[0]
        file_path = os.path.join(folder_name, file_name)

        if self.flag_cleaning_geo:
            if os.path.exists(file_path):
                os.remove(file_path)

        return file_path
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to define section tr
    def organize_section_tr(self, section_info, domain_name):

        log_stream.info(' -----> Organize section tr ... ')

        template_tags = self.alg_template_tags

        for section_tag_ref, section_fields_ref in section_info.items():

            section_name_ref = section_fields_ref['section_description']
            if 'alias' in list(section_fields_ref.keys()):
                section_alias_ref = section_fields_ref['alias']
            else:
                section_alias_ref = []

            section_name_generic = [section_name_ref]
            if section_alias_ref is not None:
                section_name_generic.extend(section_alias_ref)

            log_stream.info(' ------> Section "' + section_name_ref + '" ... ')

            if self.method_tr == 'method_regional':

                section_par_a, section_par_b, section_par_correction_factor, \
                    section_par_tr, section_par_qr = get_tr_params(section_name_ref)

                section_info[section_tag_ref]['section_par_a'] = section_par_a
                section_info[section_tag_ref]['section_par_b'] = section_par_b
                section_info[section_tag_ref]['section_par_correction_factor'] = section_par_correction_factor
                section_info[section_tag_ref]['section_par_tr'] = section_par_tr    # FIX: tr and qr are switched
                section_info[section_tag_ref]['section_par_qr'] = section_par_qr    # FIX: tr and qr are switched

            elif self.method_tr == 'method_q_t':

                file_type = self.file_type_hydro_qt
                folder_name_raw, file_name_raw = self.folder_name_hydro_qt, self.file_name_hydro_qt

                file_path_def = None
                for section_name_tmp in section_name_generic:

                    template_values = {'section_name': section_name_tmp, 'domain_name': domain_name}
                    folder_name_def = fill_tags2string(folder_name_raw, template_tags, template_values)
                    if isinstance(folder_name_def, tuple):
                        folder_name_def = folder_name_def[0]
                    file_name_def = fill_tags2string(file_name_raw, template_tags, template_values)
                    if isinstance(file_name_def, tuple):
                        file_name_def = file_name_def[0]

                    if os.path.exists(os.path.join(folder_name_def, file_name_def)):
                        file_path_def = os.path.join(folder_name_def, file_name_def)
                        break

                if file_path_def is not None:
                    if file_type is not None:
                        if file_type == 'csv':
                            section_tr_series = get_tr_file(file_path_def)
                            section_info[section_tag_ref]['section_q_t'] = section_tr_series
                        else:
                            log_stream.error(' ===> File type of tr "' + file_type + '" is not supported')
                            raise NotImplemented('Case not implemented yet')
                    else:
                        log_stream.error(' ===> File type of tr is defined by NoneType')
                        raise NotImplemented('Case not implemented yet')
                else:
                    log_stream.error(' ===> File tr was not defined by the procedure. '
                                     'Probably the file does not exist. Check the section name(s).')
                    raise RuntimeError("File tr must be defined to correctly run the algorithm")

            else:
                log_stream.error(' ===> Tr method "' + self.method_tr +
                                 '" is not available. Supported methods are: "method_regional" or "method_q_t"')
                raise NotImplemented('Case not implemented yet')

            log_stream.info(' ------> Section "' + section_name_ref + '" ... DONE')

        log_stream.info(' -----> Organize section tr ... DONE')

        return section_info

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to define upstream and downstream link(s)
    @staticmethod
    def organize_section_links(section_info):

        log_stream.info(' -----> Organize section links ... ')

        section_drainage_area_tags = ['section_drainage_area', 'drainage_area']

        for section_tag_ref, section_fields_ref in section_info.items():
            section_name_upstream_ref = section_fields_ref['name_point_upstream']
            section_name_downstream_ref = section_fields_ref['name_point_downstream']
            section_name_ref = section_fields_ref['section_description']

            log_stream.info(' ------> Section "' + section_name_ref + '" ... ')

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

                        section_field_list = join_name_and_alias(
                            section_fields_tmp, field_description='section_description', field_alias='alias')
                        section_name_tmp = section_name_step.lower()

                        for section_field_step in section_field_list:
                            section_field_step = section_field_step.lower()
                            if section_field_step == section_name_tmp:
                                section_fields_step = section_fields_tmp.copy()
                                break

                        if section_fields_step is not None:
                            break

                    if section_fields_step is None:
                        log_stream.error(' ===> Section downstream "' + section_name_step +
                                         '" is not available in the section list')
                        raise RuntimeError('Check your configuration file and remove section if not declared or needed')

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

                        section_field_list = join_name_and_alias(
                            section_fields_tmp, field_description='section_description', field_alias='alias')
                        section_name_tmp = section_name_step.lower()

                        for section_field_step in section_field_list:
                            section_field_step = section_field_step.lower()
                            if section_field_step == section_name_tmp:
                                section_fields_step = section_fields_tmp.copy()
                                break

                        if section_fields_step is not None:
                            break

                    if section_fields_step is None:
                        log_stream.error(' ===> Section upstream "' + section_name_step +
                                         '" is not available in the section list')
                        raise RuntimeError('Check your configuration file and remove section if not declared or needed')

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

            log_stream.info(' ------> Section "' + section_name_ref + '" ... DONE')

        log_stream.info(' -----> Organize section links ... DONE')

        return section_info

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to get section geo information
    @staticmethod
    def organize_section_geo(geo_data, geo_fields=None, section_key='section_{:}'):

        log_stream.info(' -----> Organize section data ... ')

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

        log_stream.info(' -----> Organize section data ... DONE')

        return section_collections

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to get section hydraulic information
    @staticmethod
    def organize_section_hydraulic(section_data, hydraulic_data, hydraulic_key_expected=None):

        log_stream.info(' -----> Organize section hydraulic ... ')

        if hydraulic_key_expected is None:
            hydraulic_key_expected = ['name_point_outlet', 'name_point_downstream',
                                      'name_point_upstream', 'name_point_obs', 'alias']

        section_data_tmp = deepcopy(section_data)
        for section_key, section_fields in section_data_tmp.items():
            section_description = section_fields['section_description']

            log_stream.info(' ------> Section "' + section_description + '" ... ')

            hydraulic_ws = None
            for hydraulic_id, hydraulic_fields in hydraulic_data.items():
                hydraulic_description = hydraulic_fields['description']

                if not isinstance(hydraulic_description, list):
                    hydraulic_description = [hydraulic_description]

                hydraulic_alias = None
                if "alias" in list(hydraulic_fields.keys()):
                    hydraulic_alias = hydraulic_fields['alias']

                if hydraulic_alias is not None:
                    hydraulic_description.extend(hydraulic_alias)

                for hydraulic_tmp in hydraulic_description:
                    if hydraulic_tmp == section_description:
                        hydraulic_ws = deepcopy(hydraulic_fields)
                        break

                if hydraulic_ws is not None:
                    break

            if hydraulic_ws is not None:
                for hydraulic_key, hydraulic_value in hydraulic_ws.items():
                    if hydraulic_key not in list(section_fields.keys()):
                        if hydraulic_key in hydraulic_key_expected:
                            section_fields[hydraulic_key] = hydraulic_value
                        if (hydraulic_key == 'drainage_area') and (section_fields['section_drainage_area'] == -9999):
                            section_fields['section_drainage_area'] = hydraulic_value

                section_data[section_key] = section_fields

                log_stream.info(' ------> Section "' + section_description + '" ... DONE')

            else:

                log_stream.info(' ------> Section "' + section_description + '" ... SKIPPED')
                log_stream.warning(' ===> Hydraulic info for section "' + section_description + '" are not available')
                section_data.pop(section_key)

        if section_data.__len__() != section_data_tmp.__len__():
            log_stream.info(' -----> Organize section hydraulic ... FAILED')
            log_stream.error(' ===> Hydraulic information are not available for all sections')
            raise RuntimeError('Hydraulic information must be available for all sections defined in the geo file')
        else:
            log_stream.info(' -----> Organize section hydraulic ... DONE')

        return section_data

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to get map geo information
    @staticmethod
    def organize_map_geo(geo_data, geo_fields=None):

        log_stream.info(' -----> Organize map data ... ')

        if geo_fields is None:
            geo_fields = ['area_reference_id', 'area_mask', 'area_mask_extended',
                          'area_reference_geo_x', 'area_reference_geo_y',
                          'coord_bottom_y', 'coord_top_y', 'coord_left_x', 'coord_right_x',
                          'area_epsg_code',
                          'area_reference_shape', 'idx_reference_unique',
                          'idx_reference_file', 'idx_reference_n']   # ADD fields to memory saver

        map_collections = {}
        for geo_field in geo_fields:

            if geo_field in list(geo_data.keys()):
                geo_value = geo_data[geo_field]

                if hasattr(geo_value, 'ndim'):
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
                        log_stream.error(' ===> Geo field "' + geo_field + '" format for array is not supported')
                        raise NotImplementedError('Case not implemented yet')
                else:
                    if isinstance(geo_value, int):
                        geo_value = str(geo_value)
                    elif isinstance(geo_value, list):
                        pass
                    elif isinstance(geo_value, tuple):
                        pass
                    else:
                        log_stream.error(' ===> Geo field "' + geo_field + '" format for scalar is not supported')
                        raise NotImplementedError('Case not implemented yet')

                map_collections[geo_field] = geo_value
            else:

                if geo_field == 'area_reference_geo_x':
                    log_stream.warning(' ===> Geo field "' + geo_field + '" will be defined by default.')
                elif geo_field == 'area_reference_geo_y':
                    log_stream.warning(' ===> Geo field "' + geo_field + '" will be defined by default.')
                elif geo_field == 'area_reference_id':
                    log_stream.warning(' ===> Geo field "' + geo_field +
                                       '" is not defined. If "area_reference_geo_x" and "area_reference_geo_y" '
                                       'are not defined, the algorithm will exit for this reason.')
                elif geo_field == 'area_reference_shape':
                    log_stream.warning(' ===> Geo field "' + geo_field +
                                       '" must be defined for applications to avoid the ram memory issues.'
                                       'Please, compute again the geo file if the field is not available.')
                elif geo_field == 'idx_reference_unique':
                    pass
                elif geo_field == 'idx_reference_file':
                    pass
                elif geo_field == 'idx_reference_n':
                    pass
                else:
                    log_stream.error(' ===> Geo field "' + geo_field + '" must be defined')
                    raise IOError('Field not found')

        if ('area_reference_geo_x' not in list(map_collections.keys())) or \
                ('area_reference_geo_y' not in list(map_collections.keys())):

            if 'area_reference_id' in list(map_collections.keys()):

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

            else:
                log_stream.error(' ===> Geo field "area_reference_id" must be defined to define'
                                 'by default "area_reference_geo_x" and "area_reference_geo_y"')
                raise RuntimeError('Field not found. Check the availability of the field "area_reference_id"')

        if 'area_epsg_code' in list(map_collections.keys()):
            epsg_code_unchecked = map_collections['area_epsg_code']
            epsg_code_checked = check_epsg_code(epsg_code_unchecked)
            map_collections['area_epsg_code'] = epsg_code_checked

        else:
            log_stream.error(' ===> Geo field "area_epsg_code" must be defined')
            raise IOError('Field not found')

        log_stream.info(' -----> Organize map data ... DONE')

        return map_collections

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to organize geographical data
    def organize_geo(self):

        # info start organize datasets
        log_stream.info(' ---> Organize geographical datasets ... ')

        # iterate over domain list
        geo_collection = {}
        for domain_name in self.domain_name_list:

            # info start domain
            log_stream.info(' ----> Domain "' + domain_name + '" ... ')

            # define file path(s)
            file_path_geo_data_collections = self.define_domain_collection(
                domain_name, self.folder_name_geo_data_out, self.file_name_geo_data_out)

            file_path_geo_area_collections = self.define_domain_collection(
                domain_name, self.folder_name_geo_area_out, self.file_name_geo_area_out)

            file_path_geo_memory_area_collections = self.define_domain_collection(
                domain_name, self.folder_name_geo_memory_area_out, self.file_name_geo_memory_area_out)

            # check memory saver option
            if self.memory_saver:
                if not os.path.exists(file_path_geo_memory_area_collections):
                    if os.path.exists(file_path_geo_data_collections):
                        os.remove(file_path_geo_data_collections)
                        os.remove(file_path_geo_area_collections)

            # check datasets availability
            if not os.path.exists(file_path_geo_data_collections):

                # get static information
                geo_data = self.read_geo_info_generic(domain_name)
                hydraulic_data = self.read_geo_info_hydraulic(domain_name)
                hydro_data = self.read_hydro_info(domain_name)

                # adjust section data from raw format
                section_data = self.organize_section_geo(geo_data)

                # remap description names (to match section and hydraulic objects)
                section_data = map_description_section_2_hydraulic(section_data, hydraulic_data)
                # remap description names (to match hydrological and hydraulic objects)
                hydro_data = map_description_hydro_2_hydraulic(hydro_data, hydraulic_data)

                # organize section datasets
                section_data = self.organize_section_hydraulic(section_data, hydraulic_data)
                section_data = self.organize_section_links(section_data)
                section_data = self.organize_section_tr(section_data, domain_name)
                # organize map datasets
                map_data = self.organize_map_geo(geo_data)

                area_data = map_data['area_reference_id']
                coord_by, coord_ty = map_data['coord_bottom_y'], map_data['coord_top_y']
                coord_rx, coord_lx = map_data['coord_right_x'], map_data['coord_left_x']
                epsg_code = map_data['area_epsg_code']

                # memory saver (organize datasets)
                if self.memory_saver:

                    var_area = deepcopy(map_data['area_reference_id'])
                    # var_idx_list = deepcopy(map_data['idx_reference_unique'])
                    map_data['area_reference_id'] = file_path_geo_memory_area_collections

                    var_idx_file, var_idx_n = {}, {}
                    for section_id, section_fields in section_data.items():
                        section_description = section_fields['section_description']
                        section_area = int(section_fields['section_area'])

                        if section_area > 0:
                            file_path_geo_memory_idx_collections = self.define_domain_collection(
                                domain_name, self.folder_name_geo_memory_idx_out, self.file_name_geo_memory_idx_out,
                                section_description=section_description)
                            var_idx_file[section_description] = file_path_geo_memory_idx_collections
                            var_idx_n[section_description] = section_area

                    map_data['idx_reference_file'] = var_idx_file
                    map_data['idx_reference_n'] = var_idx_n
                else:
                    var_area, var_idx_list, var_idx_file = None, None, None

                # store obj data
                obj_data = {'section_data': section_data, 'map_data': map_data, 'hydro_data': hydro_data}

                # create folder
                folder_name_collections, file_name_collections = os.path.split(file_path_geo_data_collections)
                os.makedirs(folder_name_collections, exist_ok=True)
                # write geographical collections
                write_obj(file_path_geo_data_collections, obj_data)

                # create folder
                folder_name_collections, file_name_collections = os.path.split(file_path_geo_area_collections)
                os.makedirs(folder_name_collections, exist_ok=True)

                write_geo_area(file_name=file_path_geo_area_collections,
                               file_data=area_data,
                               file_geo_x_west=coord_lx, file_geo_x_east=coord_rx,
                               file_geo_y_south=coord_by, file_geo_y_north=coord_ty,
                               file_geo_x=None,  # obj_data['map_data']['area_reference_geo_x'],
                               file_geo_y=None,  # obj_data['map_data']['area_reference_geo_y'],
                               file_epsg_code=epsg_code, file_metadata=None)

                # memory saver (dump datasets)
                if self.memory_saver:

                    # create folder
                    folder_name_collections, file_name_collections = os.path.split(file_path_geo_memory_area_collections)
                    os.makedirs(folder_name_collections, exist_ok=True)
                    # write ara large variable
                    write_obj(file_path_geo_memory_area_collections, var_area)

                    # iterate over idx list
                    for section_id, section_fields in section_data.items():
                        section_description = section_fields['section_description']
                        section_area = int(section_fields['section_area'])

                        # skip the first index (not defined
                        if section_area > 0:

                            # get idx file path
                            file_path_geo_memory_idx_step = var_idx_file[section_description]

                            # check data availability
                            if not os.path.exists(file_path_geo_memory_idx_step):

                                # compute idx variable
                                var_idx_data = np.argwhere(var_area == section_area)

                                # create folder
                                folder_name_collections, file_name_collections = os.path.split(
                                    file_path_geo_memory_idx_step)
                                os.makedirs(folder_name_collections, exist_ok=True)
                                # write idx large variable
                                write_obj(file_path_geo_memory_idx_step, var_idx_data)

                # info end domain
                log_stream.info(' ----> Domain "' + domain_name + '" ... DONE')

            else:

                # info end domain (data loaded from previous application)
                log_stream.info(' ----> Domain "' + domain_name + '" ... SKIPPED. Datasets previously computed.')
                obj_data = read_obj(file_path_geo_data_collections)

            # store in geo collection
            geo_collection[domain_name] = obj_data

        # info end organize datasets
        log_stream.info(' ---> Organize geographical datasets ... DONE')

        return geo_collection
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
