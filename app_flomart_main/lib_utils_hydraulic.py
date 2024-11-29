"""
Library Features:

Name:          lib_utils_hydraulic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220208'
Version:       '1.0.0'
"""

#######################################################################################
# Libraries
import logging
import json

from copy import deepcopy

from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)
#######################################################################################


# -------------------------------------------------------------------------------------
# Method to get hydraulic file information
def read_file_hydraulic(file_name, tag_domain_name='domain_name', tag_domain_sections='domain_sections'):

    with open(file_name) as file_handle:
        file_data = json.load(file_handle)

    domain_name = file_data[tag_domain_name]
    domain_section = file_data[tag_domain_sections]

    return domain_section
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to map description section and hydraulic
def map_description_section_2_hydraulic(section_data_in, hydraulic_data):

    section_data_out = {}
    for section_id, section_fields in section_data_in.items():
        section_description_tmp = section_fields['section_description']

        section_description_def, hydraulic_id, hydraulic_area = None, None, None
        for hydraulic_id, hydraulic_fields in hydraulic_data.items():

            hydraulic_description = hydraulic_fields['description']
            if 'id_area' in list(hydraulic_fields.keys()):
                hydraulic_area = hydraulic_fields['id_area']
            else:
                log_stream.error(' ===> Section "id_area" is not defined. '
                                 'For some algorithm mode the key must be defined'
                                 'in the json hydraulic file')
                raise RuntimeError('Set "id_area" key for each sections in the hydraulic file')

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

        if section_description_def is None:
            log_stream.error(' ===> Section is not available in the hydraulic description or alias')
            raise RuntimeError('Section must be included. Check the section and hydraulic names')

        # update section fields according to the hydraulic id
        section_fields['section_description'] = section_description_def
        section_fields['section_area'] = hydraulic_area
        section_data_out[hydraulic_id] = section_fields

    return section_data_out
# -------------------------------------------------------------------------------------
