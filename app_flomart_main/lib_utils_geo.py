"""
Library Features:

Name:          lib_utils_geo
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20241129'
Version:       '1.1.0'
"""

#######################################################################################
# Libraries
import logging
import re

import numpy as np

from lib_utils_io import read_mat
from lib_info_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# File Variable(s) - Version 1
# AreeCompetenza, EPSG_domain, LonArea, LatArea, LatLL, LatUR, LonLL, LonUR,
# bacino_sezione_sort, nomi_sezioni_sort, indici_sort, a1dQindex,
# mappa_aree, mappa_aree_allargata, drainage_area_section_km2

# File Variable(s) - Version 2
# AreeCompetenza, EPSG_domain, {LonArea}, {LatArea}, LatLL, LatUR, LonLL, LonUR,
# bacino_sezione_sort, nomi_sezioni_sort, indici, a1dQindex_sort,
# mappa_aree, mappa_aree_allargata, drainage_area_section_km2

# Geographical lookup table
geo_lookup_table = {'AreeCompetenza': 'area_reference_id',
                    'EPSG_domain': 'area_epsg_code',
                    'Lon_dominio_UTM': 'area_reference_geo_x', 'Lat_dominio_UTM': 'area_reference_geo_y',
                    'LatLL': 'coord_bottom_y', 'LatUR': 'coord_top_y',
                    'LonLL': 'coord_left_x', 'LonUR': 'coord_right_x',
                    'bacino_sezione_sort': 'section_description', 'nomi_sezioni_sort': 'section_name',
                    'indici': 'section_id', 'a1dQindex_sort': 'section_discharge_idx',
                    'mappa_aree': 'area_mask', 'mappa_aree_allargata': 'area_mask_extended',
                    'drainage_area_section_km2': 'section_drainage_area'}
# Geographical type table
geo_type_table = {'AreeCompetenza': 'ancillary',
                  'EPSG_domain': 'mandatory',
                  'Lon_dominio_UTM': 'mandatory', 'Lat_dominio_UTM': 'mandatory',
                  'LatLL': 'mandatory', 'LatUR': 'mandatory',
                  'LonLL': 'mandatory', 'LonUR': 'mandatory',
                  'bacino_sezione_sort': 'mandatory', 'nomi_sezioni_sort': 'mandatory',
                  'indici': 'mandatory', 'a1dQindex_sort': 'mandatory',
                  'mappa_aree': 'mandatory', 'mappa_aree_allargata': 'mandatory',
                  'drainage_area_section_km2': 'mandatory'}
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read geographical file in mat format
def read_file_geo(file_name, excluded_fields=None):

    if excluded_fields is None:
        excluded_fields = ['__header__', '__version__', '__globals__']

    if file_name.endswith('mat'):

        # read geo datasets
        file_data = read_mat(file_name)

        # check the keys in type and the lut dict(s)
        key_list_lut = sorted(list(geo_lookup_table.keys()))
        key_list_type = sorted(list(geo_type_table.keys()))
        if key_list_lut != key_list_type:
            log_stream.error(' ===> Geographical LUT and Type lists are not the same')
            raise RuntimeError('The keys must be the same in both lists')

        # check fields
        for file_key, file_type in geo_type_table.items():
            if file_key not in list(file_data.keys()):
                if file_type == 'mandatory':
                    log_stream.error(' ===> File mandatory field "' + file_key +
                                     '" is not find in the geographical reference file "' + file_name + '"')
                    raise RuntimeError('Field is mandatory and must be in the geographical reference file')
                elif file_type == 'ancillary':
                    log_stream.warning(' ===> File ancillary field "' + file_key +
                                       '" is not find in the geographical reference file "' + file_name + '"')
                else:
                    log_stream.error(' ===> File type "' + file_type + '" is not supported')
                    raise NotImplementedError('Case not implemented yet')

        # organize fields
        file_collections = {}
        for file_key, file_values in file_data.items():

            if (file_key not in excluded_fields) and (file_key in list(geo_lookup_table.keys())):
                var_idx = list(geo_lookup_table.keys()).index(file_key)
                var_name = list(geo_lookup_table.values())[var_idx]
                file_collections[var_name] = file_values

                # to mitigate the amount of memory for "area_reference_id" field
                if var_name == 'area_reference_id':
                    file_unique = np.unique(file_values)
                    file_unique = list(file_unique.astype(int))
                    file_shape = file_values.shape
                    file_collections['idx_reference_unique'] = file_unique
                    file_collections['area_reference_shape'] = file_shape

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
        elif isinstance(epsg_parts, list) and (epsg_parts.__len__() == 2):
            epsg_value = epsg_parts[0]
        else:
            log_stream.error(' ===> EPSG format is not supported')
            raise NotImplementedError('Case not implemented yet')

        epsg_code = epsg_prefix + epsg_sep + epsg_value

    return epsg_code

# -------------------------------------------------------------------------------------
