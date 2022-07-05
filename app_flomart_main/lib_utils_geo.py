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
import numpy as np


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
                    'LonArea': 'area_reference_geo_x', 'LatArea': 'area_reference_geo_y',
                    'LatLL': 'coord_bottom_y', 'LatUR': 'coord_top_y',
                    'LonLL': 'coord_left_x', 'LonUR': 'coord_right_x',
                    'bacino_sezione_sort': 'section_description', 'nomi_sezioni_sort': 'section_name',
                    'indici_sort': 'section_id', 'a1dQindex': 'section_discharge_idx',
                    'mappa_aree': 'area_mask', 'mappa_aree_allargata': 'area_mask_extended'}
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
# Method to convert decimal degrees to km [method2]
def deg_2_km_OLD(deg, lat=None):
    if lat is None:
        km = deg * 110.54
    else:
        km = deg * 111.32 * np.cos(np.deg2rad(lat))
    return km
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to compute cell area in m^2
def compute_cell_area_OLD(geo_x, geo_y, cell_size_x, cell_size_y):

    # Method constant(s)
    r = 6378388  # (Radius)
    e = 0.00672267  # (Ellipsoid)

    # dx = (R * cos(lat)) / (sqrt(1 - e2 * sqr(sin(lat)))) * PI / 180
    dx_2d = (r * np.cos(np.abs(geo_y) * np.pi / 180)) / \
            (np.sqrt(1 - e * np.sqrt(np.sin(np.abs(geo_y) * np.pi / 180)))) * np.pi / 180
    # dy = (R * (1 - e2)) / pow((1 - e2 * sqr(sin(lat))),1.5) * PI / 180
    dy_2d = (r * (1 - e)) / np.power((1 - e * np.sqrt(np.sin(np.abs(geo_y) / 180))), 1.5) * np.pi / 180

    # area cell in m^2
    area_cell = ((dx_2d / (1 / cell_size_x)) * (dy_2d / (1 / cell_size_y)))  # [m^2]

    return area_cell

# -------------------------------------------------------------------------------------
