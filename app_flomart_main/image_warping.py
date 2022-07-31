import os
from osgeo import gdal
import rasterio

files_to_mosaic = [
    "/home/fabio/Desktop/PyCharm_Workspace/flomart-ws/data_static/hazard_data/Abaco_Caccamo-Grazie/Caccamo-Grazie_WD_max_T020.tif",
    "/home/fabio/Desktop/PyCharm_Workspace/flomart-ws/data_static/hazard_data/Abaco_Grazie-Foce/Grazie-Foce_WD_max_T020.tif",
    "/home/fabio/Desktop/PyCharm_Workspace/flomart-ws/data_static/hazard_data/Abaco_Polverina-Caccamo/Polverina-Caccamo_WD_max_T020.tif"
]

file_output = "/home/fabio/Desktop/PyCharm_Workspace/flomart-ws/data_static/hazard_data/chienti_merged_threemap_warp.tif"
file_output_2 = "/home/fabio/Desktop/PyCharm_Workspace/flomart-ws/data_static/hazard_data/chienti_merged_threemap_gdal_merge.tif"

if os.path.exists(file_output):
    os.remove(file_output)
if os.path.exists(file_output_2):
    os.remove(file_output_2)

for file_name in files_to_mosaic:
    file_handle = rasterio.open(file_name)

    file_transform = file_handle.transform
    file_meta = file_handle.meta
    print(file_meta)

g = gdal.Warp(file_output, files_to_mosaic, format="GTiff",
              options=["COMPRESS=LZW", "TILED=YES"])
g = None

files_string = " ".join(files_to_mosaic)

command = "gdal_merge.py -ps 2.0 2.0 -n 0 -o " + file_output_2 + " -of gtiff " + files_string
print(os.popen(command).read())

print('cioa')
