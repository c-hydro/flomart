import numpy as np
import matplotlib.pyplot as plt
import os, glob
from osgeo import gdal, gdalconst
import matplotlib.pyplot as plt
import subprocess

import rasterio
from rasterio.merge import merge
from rasterio.transform import Affine

for tr in range(1,501):
    print(tr)
    # sPathHazardData = '/media/matteo/Elements/LINUX/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data'
    sPathHazardData = '/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data'
    path_file_hazard_read1 = sPathHazardData + '/Polverina-Caccamo/Polverina-Caccamo_WD_max_T' + (str(tr).zfill(3)) + '.tif'
    path_file_hazard_read2 = sPathHazardData + '/Caccamo-Grazie/Caccamo-Grazie_WD_max_T' + (str(tr).zfill(3)) + '.tif'
    path_file_hazard_read3 = sPathHazardData + '/Grazie-Foce/Grazie-Foce_WD_max_T' + (str(tr).zfill(3)) +'.tif'
    output_file = sPathHazardData + '/Chienti/Chienti_not_compressed.tif'
    output_file_compressed = sPathHazardData + '/Chienti/Chienti_WD_max_T' + (str(tr).zfill(3)) + '.tif'

    if not os.path.exists(path_file_hazard_read1):
        print("ERROR: hazard map path 1 not found !!!")

    if not os.path.exists(path_file_hazard_read2):
        print("ERROR: hazard map path 2 not found !!!")

    if not os.path.exists(path_file_hazard_read3):
        print("ERROR: hazard map path 3 not found !!!")


    if os.path.exists(output_file):
        os.remove(output_file)
        print("Output file has been deleted successfully")

    if os.path.exists(output_file_compressed):
        os.remove(output_file_compressed)
        print("Output compressed file has been deleted successfully")


    files_to_mosaic = [path_file_hazard_read1 , path_file_hazard_read2, path_file_hazard_read3]
    files_string = " ".join(files_to_mosaic)
    command = "gdal_merge.py -n 0  -co COMPRESS=DEFLATE -o " + output_file_compressed + " -of gtiff " + files_string
    print(os.popen(command).read())


    # compression:
    # path_file_hazard_read_tocompress = sPathHazardData + '/Chienti/compressed/Chienti_WD_max_T' + (str(tr).zfill(3)) + '.tif'
    # gdalinput = path_file_hazard_read_tocompress
    # gdaloutput = output_file_compressed
    # translateoptions = gdal.TranslateOptions(gdal.ParseCommandLine("-of Gtiff -co COMPRESS=DEFLATE"))
    # gdal.Translate(gdaloutput, gdalinput, options=translateoptions)

print('All done')





# #################################################################
# # read:
# # file_data, file_proj, file_geotrans = read_file_tif(file_name)
# # one raster:
# file_handle = rasterio.open(files_to_mosaic[0])
# file_proj = file_handle.crs.wkt
# file_geotrans = file_handle.transform
# file_tags = file_handle.tags()
# file_bands = file_handle.count
# file_metadata = file_handle.profile
#
# if file_bands == 1:
#     file_data = file_handle.read(1)
# elif file_bands > 1:
#     file_data = []
#     for band_id in range(0, file_bands):
#         file_data_tmp = file_handle.read(band_id + 1)
#         file_data.append(file_data_tmp)
# else:
#     log_stream.error(' ===> File multi-band are not supported')
#     raise NotImplementedError('File multi-band not implemented yet')
# ########
#
# # file_data = read_file_tiff(domain_map_file_ancillary)
# plt.figure()
# plt.imshow(file_data)
# plt.colorbar()
# plt.clim(0, 8)
# plt.show()
#
#
#
#
#
#
#
#
#
# ##############
# # save tiff:
# # save_file_tiff(file_path_scenario_plot_tiff,
# #                domain_map_data, domain_geo_x, domain_geo_y,
# #                file_epsg_code=domain_epsg_code)
# # def save_file_tiff(file_name, file_data, file_geo_x, file_geo_y, file_metadata=None, file_epsg_code='EPSG:32632'):
#
# # define corners of final map:
# xmin =[]
# xmax=[]
# ymin =[]
# ymax=[]
#
# # for file in files_to_mosaic:
# #     src = gdal.Open(file)
# #     ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
# #     lrx = ulx + (src.RasterXSize * xres)
# #     lry = uly + (src.RasterYSize * yres)
# #     xmin.append(min(lrx, ulx))
# #     xmax.append(max(lrx, ulx))
# #     ymin.append(min(lry, uly))
# #     ymax.append(max(lry, uly))
#
# src = gdal.Open(files_to_mosaic[0])
# ulx, xres, xskew, uly, yskew, yres  = src.GetGeoTransform()
# lrx = ulx + (src.RasterXSize * xres)
# lry = uly + (src.RasterYSize * yres)
# xmin.append(min(lrx, ulx))
# xmax.append(max(lrx, ulx))
# ymin.append(min(lry, uly))
# ymax.append(max(lry, uly))
#
# file_data_height, file_data_width = file_data.shape
#
#
# # corners of template:
# file_geo_x_west = np.min(xmin)
# file_geo_x_east = np.max(xmax)
# file_geo_y_south = np.min(ymin)
# file_geo_y_north = np.max(ymax)
#
# file_data_transform = rasterio.transform.from_bounds(
#     file_geo_x_west, file_geo_y_south, file_geo_x_east, file_geo_y_north,
#     file_data_width, file_data_height)
#
# file_epsg_code='EPSG:32633'
#
#
#
#
#
#
#
# #########
# # write_file_tif(file_name, file_data,
# #                file_data_width, file_data_height, file_data_transform, file_epsg_code,
# #                file_metadata=file_metadata)
#
# # def write_file_tif(file_name, file_data, file_wide, file_high, file_geotrans, file_proj,
# #                    file_metadata=None,
# #                    file_format=gdalconst.GDT_Float32):
#
#
# if not isinstance(file_data, list):
#     file_data = [file_data]
#
# if file_metadata is None:
#     file_metadata = {'description_field': 'data'}
#
# if not isinstance(file_metadata, list):
#     file_metadata = [file_metadata] * file_data.__len__()
#
# if isinstance(file_data_transform, Affine):
#     file_data_transform = file_data_transform.to_gdal()
#
# file_crs = rasterio.crs.CRS.from_string(file_epsg_code)
# file_wkt = file_crs.to_wkt()
# file_n = file_data.__len__()
# file_format= gdalconst.GDT_Float32
# dset_handle = gdal.GetDriverByName('GTiff').Create("/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data/merged2.tif",
#                                                    file_data_width, file_data_height, file_n, file_format,
#                                                    options=['COMPRESS=DEFLATE'])
# dset_handle.SetGeoTransform(file_data_transform)
# dset_handle.SetProjection(file_wkt)
#
#
#
# # for file_id, (file_data_step, file_metadata_step) in enumerate(zip(file_data, file_metadata)):
# #     print(file_id)
# #     dset_handle.GetRasterBand(file_id + 1).WriteArray(file_data_step)
# #     dset_handle.GetRasterBand(file_id + 1).SetMetadata(file_metadata_step)
#
# dset_handle.GetRasterBand(1).WriteArray(file_data[0])
# #dset_handle.GetRasterBand(1).SetMetadata(file_metadata[0])
#
# del dset_handle
# # -------------------
#
#
#
#
#
# # list = ["/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data/merged2.tif",
# #         "/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data/merged2.tif",
# #        "/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data/merged2.tif"]
# # raster_to_mosaic =[]
# # for p in list:
# #     raster = rasterio.open(p)
# #     raster_to_mosaic.append(raster)
#
#
#
#
#
#
#
#
# command = "gdal_merge.py -n 0 -o /home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data/mergedmap.tif -of gtiff " + files_string
# print(os.popen(command).read())
#
#
# from osgeo import gdal
# g = gdal.Warp("/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data/mergedmap.tif",
#               files_to_mosaic, format="GTiff",options=["COMPRESS=LZW", "TILED=YES"] )
#
#
#
# with rasterio.open(files_to_mosaic[0]) as src:
#     meta = src.meta.copy()
#
# # The merge function returns a single array and the affine transform info
# arr, out_trans = merge(files_to_mosaic)
#
# meta.update({
#     "driver": "GTiff",
#     "height": arr.shape[1],
#     "width": arr.shape[2],
#     "transform": out_trans
# })
#
# # Write the mosaic raster to disk
# output = "/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data/mergedmap.tif"
# with rasterio.open(output, "w", **meta) as dest:
#     dest.write(arr)
#
#
#
# cmd = 'gdal_merge.py -ps 10 -10 -o /home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data/mergedmap.tif'
# subprocess.call(cmd.split() + files_to_mosaic)
#
# vrt = gdal.BuildVRT("/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data/merged.vrt", files_to_mosaic)
# gdal.Translate("/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data/merged2.tif", vrt)
# vrt =None
#
#
# noNaN_path = "/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data/"
# ndvi_name = "merged"
#
# for pix in files_to_mosaic:
#     gdal_cmd = 'gdal_translate -a_nodata nan {0} {1}/{2}_{3}.tif'.format(pix, noNaN_path, ndvi_name, files_to_mosaic.index(pix))
#     subprocess.call(gdal_cmd, shell=True)
#
# subprocess.call(gdal_cmd, shell=True)
#
# vrt = gdal.BuildVRT("/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data/merged.vrt", files_to_mosaic,
#                     VRTNodata="nan", srcNodata="nan")
# gdal.Translate("/home/matteo/Documents/CIMA_projects/RT_FloodMapping/data/data_Marche_Chienti/data_static/hazard_data/merged2.tif", vrt)
#
#
