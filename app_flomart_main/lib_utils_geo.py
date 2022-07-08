"""
Library Features:

Name:          lib_utils_geo
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220208'
Version:       '1.0.0'
"""

#######################################################################################
# Libraries
import logging
import re

from lib_utils_io import read_mat
from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Geographical lookup table
geo_lookup_table = {'AreeCompetenza': 'area_reference_id',
                    'EPSG_domain': 'area_epsg_code',
                    'LonArea': 'area_reference_geo_x', 'LatArea': 'area_reference_geo_y',
                    'LatLL': 'coord_bottom_y', 'LatUR': 'coord_top_y',
                    'LonLL': 'coord_left_x', 'LonUR': 'coord_right_x',
                    'bacino_sezione_sort': 'section_description', 'nomi_sezioni_sort': 'section_name',
                    'indici_sort': 'section_id', 'a1dQindex': 'section_discharge_idx',
                    'mappa_aree': 'area_mask', 'mappa_aree_allargata': 'area_mask_extended',
                    'drainage_area_section_km2': 'section_drainage_area'}
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read geographical file in mat format
def read_file_geo(file_name, excluded_fields=None):

    if excluded_fields is None:
        excluded_fields = ['__header__', '__version__', '__globals__']

    if file_name.endswith('mat'):

        file_data = read_mat(file_name)

        file_collections = {}
        for file_key, file_values in file_data.items():
            if (file_key not in excluded_fields) and (file_key in list(geo_lookup_table.keys())):
                var_idx = list(geo_lookup_table.keys()).index(file_key)
                var_name = list(geo_lookup_table.values())[var_idx]
                file_collections[var_name] = file_values
            elif (file_key in excluded_fields) and (file_key not in list(geo_lookup_table.keys())):
                pass
            elif (file_key not in excluded_fields) and (file_key not in list(geo_lookup_table.keys())):
                pass
            else:
                log_stream.error(' ===> Field expected "' + file_key + '" is not found')
                raise IOError('Field must be set in the "' + file_name + '"')
    else:
        log_stream.error(' ===> File format "' + file_name + '" is not supported')
        raise NotImplementedError('Case not implemented yet')

    return file_collections
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to check epsg code
def check_epsg_code(epsg_value, epsg_prefix='EPSG', epsg_sep=':'):

    if not isinstance(epsg_value, str):
        epsg_value = str(epsg_value)

    if epsg_value.isnumeric():
        epsg_code = epsg_prefix + epsg_sep + epsg_value
    else:
        epsg_parts = re.findall(r'\d+', epsg_value)

        if isinstance(epsg_parts, list) and (epsg_parts.__len__() == 1):
            epsg_value = epsg_parts[0]
        else:
            log_stream.error(' ===> EPSG format is not supported')
            raise NotImplementedError('Case not implemented yet')

        epsg_code = epsg_prefix + epsg_sep + epsg_value

    return epsg_code

# -------------------------------------------------------------------------------------
