"""
Library Features:

Name:          lib_utils_io
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220208'
Version:       '1.0.0'
"""

#######################################################################################
# Libraries
import logging
import warnings
import tempfile
import os
import json
import pickle
import rasterio
import numpy as np
import xarray as xr
import pandas as pd
import scipy.io
try:
    import h5py
except ImportError:
    warnings.warn(" ===> H5py library is not imported. File .mat will not correctly read")

from copy import deepcopy
from rasterio.transform import Affine
from osgeo import gdal, gdalconst
from lib_info_args import logger_name, time_format_algorithm

logging.getLogger('rasterio').setLevel(logging.WARNING)

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
import matplotlib.pylab as plt
#######################################################################################


# -------------------------------------------------------------------------------------
# Method to create a tmp name
def create_filename_tmp(prefix='tmp_', suffix='.tiff', folder=None):

    if folder is None:
        folder = '/tmp'

    with tempfile.NamedTemporaryFile(dir=folder, prefix=prefix, suffix=suffix, delete=False) as tmp:
        temp_file_name = tmp.name
    return temp_file_name
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to write file tiff
def write_file_tif(file_name, file_data, file_wide, file_high, file_geotrans, file_proj,
                   file_metadata=None,
                   file_format=gdalconst.GDT_Float32):

    if not isinstance(file_data, list):
        file_data = [file_data]

    if file_metadata is None:
        file_metadata = {'description_field': 'data'}
    if not isinstance(file_metadata, list):
        file_metadata = [file_metadata] * file_data.__len__()

    if isinstance(file_geotrans, Affine):
        file_geotrans = file_geotrans.to_gdal()

    if isinstance(file_proj, str):
        file_crs = rasterio.crs.CRS.from_string(file_proj)
    elif isinstance(file_proj, int):
        file_crs = rasterio.crs.CRS.from_epsg(file_proj)
    else:
        log_stream.error(' ===> Projection format is not allowed')
        raise NotImplementedError('Case not implemented yet')
    file_wkt = file_crs.to_wkt()

    file_n = file_data.__len__()
    dset_handle = gdal.GetDriverByName('GTiff').Create(file_name, file_wide, file_high, file_n, file_format,
                                                       options=['COMPRESS=DEFLATE'])
    dset_handle.SetGeoTransform(file_geotrans)
    dset_handle.SetProjection(file_wkt)

    for file_id, (file_data_step, file_metadata_step) in enumerate(zip(file_data, file_metadata)):
        dset_handle.GetRasterBand(file_id + 1).WriteArray(file_data_step)
        dset_handle.GetRasterBand(file_id + 1).SetMetadata(file_metadata_step)
    del dset_handle
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read data values in geotiff format
def read_file_tif(file_name):
    file_handle = rasterio.open(file_name)
    file_proj = file_handle.crs.wkt
    file_geotrans = file_handle.transform

    file_tags = file_handle.tags()
    file_bands = file_handle.count
    file_metadata = file_handle.profile

    if file_bands == 1:
        file_data = file_handle.read(1)
    elif file_bands > 1:
        file_data = []
        for band_id in range(0, file_bands):
            file_data_tmp = file_handle.read(band_id + 1)
            file_data.append(file_data_tmp)
    else:
        log_stream.error(' ===> File multi-band are not supported')
        raise NotImplementedError('File multi-band not implemented yet')

    return file_data, file_proj, file_geotrans
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to save data info in json format
def save_file_json(file_name, file_data_dict, file_indent=4, file_sep=',', file_float_decimals=2, file_nodata=-9999.0):

    file_data_json = {}
    for file_key, file_value in file_data_dict.items():

        if isinstance(file_value, list):
            file_value = [str(i) for i in file_value]
            file_value = file_sep.join(file_value)
        elif isinstance(file_value, int):
            if np.isnan(value_step):
                value_step = file_nodata
            file_value = str(int(file_value))
        elif isinstance(file_value, float):
            if np.isnan(value_step):
                value_step = file_nodata
            file_value = str(round(file_value, file_float_decimals))
        elif isinstance(file_value, str):
            pass
        elif isinstance(file_value, dict):
            file_tmp = {}
            for value_key, value_data in file_value.items():
                if isinstance(value_data, np.datetime64):
                    time_stamp = pd.to_datetime(str(value_data))
                    time_str = time_stamp.strftime(time_format_algorithm)
                    file_tmp[value_key] = time_str
                elif isinstance(value_data, list):
                    list_obj = []
                    for value_step in value_data:

                        if isinstance(value_step, pd.Timestamp):
                            value_obj = value_step.strftime(time_format_algorithm)
                        elif isinstance(value_step, int):
                            if np.isnan(value_step):
                                value_step = file_nodata
                            value_obj = str(int(value_step))
                        elif isinstance(value_step, float):
                            if np.isnan(value_step):
                                value_step = file_nodata
                            value_obj = str(round(value_step, file_float_decimals))
                        elif isinstance(value_step, str):
                            value_obj = deepcopy(value_step)
                        else:
                            log_stream.error(' ===> Type of list element is not supported')
                            raise NotImplementedError('Format not implemented yet')
                        list_obj.append(value_obj)
                    string_obj = file_sep.join(list_obj)
                    file_tmp[value_key] = string_obj
                elif isinstance(value_data, dict):

                    for key_step, value_step in value_data.items():
                        if isinstance(value_step, pd.Timestamp):
                            value_obj = value_step.strftime(time_format_algorithm)
                        elif isinstance(value_step, bool):
                            value_obj = str(deepcopy(value_step))
                        elif isinstance(value_step, int):
                            if np.isnan(value_step):
                                value_step = file_nodata
                            value_obj = str(int(value_step))
                        elif isinstance(value_step, float):
                            if np.isnan(value_step):
                                value_step = file_nodata
                            value_obj = str(round(value_step, file_float_decimals))
                        elif isinstance(value_step, str):
                            value_obj = deepcopy(value_step)
                        elif isinstance(value_step, list):
                            value_tmp = [str(i) for i in value_step]
                            value_obj = file_sep.join(value_tmp)
                        else:
                            log_stream.error(' ===> Type of dict element is not supported')
                            raise NotImplementedError('Format not implemented yet')

                        file_tmp[key_step] = value_obj
                else:
                    file_tmp[value_key] = value_data
            file_value = deepcopy(file_tmp)
        elif file_value is None:
            log_stream.warning(' ===> Datasets obj is defined by NoneType. Check if datasets are available or not')
        else:
            log_stream.error(' ===> Error in getting datasets')
            raise RuntimeError('Datasets case not implemented yet')

        file_data_json[file_key] = file_value

    file_data = json.dumps(file_data_json, indent=file_indent, ensure_ascii=False, sort_keys=True)
    with open(file_name, "w", encoding='utf-8') as file_handle:
        file_handle.write(file_data)

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read file json
def read_file_json(file_name):
    with open(file_name) as file_handle:
        file_data = json.load(file_handle)  #, encoding='utf-8')
    return file_data
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to create a data array
def create_darray_3d(data, time, geo_x, geo_y, geo_1d=True,
                     coord_name_x='west_east', coord_name_y='south_north', coord_name_time='time',
                     dim_name_x='west_east', dim_name_y='south_north', dim_name_time='time',
                     dims_order=None):

    if dims_order is None:
        dims_order = [dim_name_y, dim_name_x, dim_name_time]

    if geo_1d:
        if geo_x.shape.__len__() == 2:
            geo_x = geo_x[0, :]
        if geo_y.shape.__len__() == 2:
            geo_y = geo_y[:, 0]

        data_da = xr.DataArray(data,
                               dims=dims_order,
                               coords={coord_name_time: (dim_name_time, time),
                                       coord_name_x: (dim_name_x, geo_x),
                                       coord_name_y: (dim_name_y, geo_y)})
    else:
        log_stream.error(' ===> Longitude and Latitude must be 1d')
        raise IOError('Variable shape is not valid')

    return data_da
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to create a data array
def create_darray_2d(data, geo_x, geo_y, geo_1d=True, name='geo',
                     coord_name_x='west_east', coord_name_y='south_north',
                     dim_name_x='west_east', dim_name_y='south_north',
                     dims_order=None):

    if dims_order is None:
        dims_order = [dim_name_y, dim_name_x]

    if geo_1d:
        if geo_x.shape.__len__() == 2:
            geo_x = geo_x[0, :]
        if geo_y.shape.__len__() == 2:
            geo_y = geo_y[:, 0]

        data_da = xr.DataArray(data,
                               dims=dims_order,
                               coords={coord_name_x: (dim_name_x, geo_x),
                                       coord_name_y: (dim_name_y, geo_y)},
                               name=name)
    else:
        log_stream.error(' ===> Longitude and Latitude must be 1d')
        raise IOError('Variable shape is not valid')

    return data_da
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read data obj
def read_obj(file_name):
    if os.path.exists(file_name):
        data = pickle.load(open(file_name, "rb"))
    else:
        data = None
    return data
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to write data obj
def write_obj(file_name, data):
    if os.path.exists(file_name):
        os.remove(file_name)
    with open(file_name, 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read mat obj
def read_mat(file_name):
    if os.path.exists(file_name):
        try:
            data = scipy.io.loadmat(file_name)

        except BaseException as base_exp:

            # mat file has been creation wit6h version 7.3 of matlab save
            log_stream.warning(
                ' ===> Geo file in MAT format "' + file_name +
                '" is not supported by scipy library with exception "' + str(base_exp) + '".')
            log_stream.warning(' ===> Try to use h5py to read it')

            file_handle = h5py.File(file_name, 'r')
            data = {}
            for key in list(file_handle.keys()):
                fields = file_handle[key]
                if isinstance(fields, h5py._hl.dataset.Dataset):
                    obj = fields[()]

                    if isinstance(obj[0][0], h5py.Reference):
                        obj_list = []
                        for obj_tmp in obj:
                            obj_name = obj_tmp[0]
                            obj_parts = file_handle[obj_name]
                            obj_value = obj_parts[()]
                            obj_dtype = obj_parts.regionref.id.dtype

                            if obj_value.shape[0] == 1 and obj_dtype == np.float64:
                                obj_tmp = obj_parts[()][0]
                                if obj_tmp.shape[0] == 1:
                                    obj_tmp = obj_tmp[0]
                                obj_list.append(obj_tmp)
                            elif obj_value.shape[0] >= 1 and obj_dtype == np.uint16:
                                obj_joins = ''.join(chr(int(i)) for i in obj_parts[:])
                                obj_list.append(obj_joins)
                            else:
                                log_stream.error(' ===> Field parser is not supported')
                                raise NotImplementedError('Case not implemented yet')
                        data[key] = obj_list
                    else:
                        if key == 'EPSG_domain':
                            obj_joins = ''.join(chr(int(i)) for i in obj[:])
                            data[key] = int(obj_joins)
                        elif key == 'indici_sort' or key == 'a1dQindex':
                            elem_list = []
                            for elem_step in obj:
                                if elem_step.shape[0] == 1:
                                    elem_step = elem_step[0]
                                elem_list.append(elem_step)
                            data[key] = elem_list
                        elif key == 'AreeCompetenza':
                            if obj.ndim == 2:
                                obj_tmp = np.transpose(obj)
                                # obj_tmp = obj.copy()
                                data[key] = obj_tmp
                            else:
                                log_stream.error(' ===> Obj must be in 2d format')
                                raise NotImplementedError('Case not implemented yet')
                        else:
                            data[key] = obj
    else:
        data = None

    return data
# -------------------------------------------------------------------------------------


