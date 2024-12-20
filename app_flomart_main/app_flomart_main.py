#!/usr/bin/python3
"""
FLOMART - FLOod MAp in Real Time

__date__ = '20241215'
__version__ = '2.2.1'
__author__ =
        'Fabio Delogu (fabio.delogu@cimafoundation.org)',
        'Matteo Darienzo (matteo.darienzo@cimafoundation.org)',
        'Flavio Pignone (flavio.pignone@cimafoundation.org)',
        'Rocco Masi (rocco.masi@cimafoundation.org)',
        'Lorenzo Campo (lorenzo.campo@cimafoundation.org)',
        'Francesco Silvestro (francesco.silvestro@cimafoundation.org)'

__library__ = 'flomart'

General command line:
python app_flomart_main.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20241215 (2.2.1) --> Add the method to dump reference areas, fix bugs in geographical destination boundaries
20241129 (2.2.0) --> Add safe memory mode to manage the ram memory usage, extend codes to manage the info domain file, fix bugs
20240530 (2.1.1) --> Fix bugs to manage the case of missing data in the observed dataset (no links mode)
20231106 (2.1.0) --> Add mode using Q_T file(s), bugs fixing and code refactoring (to adapt new sources and format)
20221215 (2.0.5) --> FLOMART refactor for liguria operational chain
20221212 (2.0.4) --> FLOMART refactor for marche operational chain
20220801 (2.0.3) --> FLOMART refactor of part to read geographical data in different format (matlab, hdf5)
20220728 (2.0.2) --> FLOMART refactor of part to preprocess domain datasets
20220715 (2.0.1) --> FLOMART refactor to part to compute different domain and settings
20220701 (2.0.0) --> FLOMART application(s) and code refactoring based on file settings and static and dynamic datasets
20220228 (1.9.1) --> Bugs fixing for the weighted scenarios part
20220223 (1.9.0) --> Bugs fixing and code updating based on pre-operational release
20220201 (1.8.0) --> Pre-operational release
20211005 (1.7.0) --> Generic release for correcting bugs and managing the empty datasets for obs/mod discharges
20210515 (1.6.0) --> Generic release for updating tools and modules
20201214 (1.5.0) --> Test release
20200522 (1.0.0) --> Beta release
"""

# -------------------------------------------------------------------------------------
# Complete library
import logging
import time
import os

from driver_data_io_geo import DriverGeo
from driver_data_io_source import DriverDischarge
from driver_data_io_destination import DriverScenario

from argparse import ArgumentParser

from lib_utils_logging import set_logging_file
from lib_utils_io import read_file_json
from lib_utils_time import set_time
from lib_info_args import logger_name, logger_format, time_format_algorithm

# Logging
log_stream = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
alg_version = '2.2.1'
alg_release = '2024-12-15'
alg_name = 'FLOMART'
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Script Main
def main():

    # -------------------------------------------------------------------------------------
    # Get algorithm settings
    alg_settings, alg_time = get_args()

    # Set algorithm settings
    data_settings = read_file_json(alg_settings)

    # Set algorithm logging
    set_logging_file(
        logger_name=logger_name,
        logger_formatter=logger_format,
        logger_file=os.path.join(data_settings['log']['folder_name'], data_settings['log']['file_name']))
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Info algorithm
    log_stream.info(' ============================================================================ ')
    log_stream.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    log_stream.info(' ==> START ... ')
    log_stream.info(' ')

    # Time algorithm information
    start_time = time.time()
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Organize time run
    time_now, time_run, time_range = set_time(
        time_run_args=alg_time,
        time_run_file=data_settings['time']['time_now'],
        time_format=time_format_algorithm,
        time_period=data_settings['time']['time_period'],
        time_frequency=data_settings['time']['time_frequency'],
        time_rounding=data_settings['time']['time_rounding']
    )
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Geographical datasets
    driver_data_geo = DriverGeo(
        src_dict=data_settings['data']['static']['source'],
        dst_dict=data_settings['data']['static']['destination'],
        alg_ancillary=data_settings['algorithm']['ancillary'],
        alg_template_tags=data_settings['algorithm']['template'],
        flag_cleaning_geo=data_settings['algorithm']['flags']['cleaning_static_data'])
    geo_data_collection = driver_data_geo.organize_geo()
    # -------------------------------------------------------------------------------------

    # ----------------------------------------------------------------------------------
    # Iterate over time range
    for time_step in time_range:

        # -------------------------------------------------------------------------------------
        # Discharge datasets
        driver_data_source_discharge = DriverDischarge(
            time_now=time_now,
            time_run=time_step,
            geo_data_collection=geo_data_collection,
            src_dict=data_settings['data']['dynamic']['source'],
            anc_dict=data_settings['data']['dynamic']['ancillary'],
            alg_ancillary=data_settings['algorithm']['ancillary'],
            alg_template_tags=data_settings['algorithm']['template'],
            flag_cleaning_anc_discharge_obs=data_settings['algorithm']['flags']['cleaning_ancillary_data_discharge_obs'],
            flag_cleaning_anc_discharge_sim=data_settings['algorithm']['flags']['cleaning_ancillary_data_discharge_sim'])
        discharge_data_collection = driver_data_source_discharge.organize_discharge()

        # Scenario datasets
        driver_data_destination_scenario = DriverScenario(
            time_now=time_now,
            time_run=time_step,
            discharge_data_collection=discharge_data_collection,
            geo_data_collection=geo_data_collection,
            src_dict=data_settings['data']['static']['source'],
            anc_dict=data_settings['data']['dynamic']['ancillary'],
            dst_dict=data_settings['data']['dynamic']['destination'],
            alg_ancillary=data_settings['algorithm']['ancillary'],
            alg_template_tags=data_settings['algorithm']['template'],
            flag_cleaning_anc_scenario_info=data_settings['algorithm']['flags']['cleaning_ancillary_data_scenario_info'],
            flag_cleaning_anc_scenario_file=data_settings['algorithm']['flags']['cleaning_ancillary_data_scenario_file'],
            flag_cleaning_anc_scenario_map=data_settings['algorithm']['flags']['cleaning_ancillary_data_scenario_maps'],
            flag_cleaning_plot_scenario=data_settings['algorithm']['flags']['cleaning_dynamic_plot'],
            flag_cleaning_data_scenario=data_settings['algorithm']['flags']['cleaning_dynamic_data']
        )
        scenario_info_collection = driver_data_destination_scenario.organize_scenario_datasets()
        scenario_map_collection, scenario_dframe_collection = driver_data_destination_scenario.compute_scenario_map(
            scenario_info_collection)
        driver_data_destination_scenario.dump_scenario_map(scenario_map_collection,
                                                           scenario_info_collection, scenario_dframe_collection)
        # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Info algorithm
    time_elapsed = round(time.time() - start_time, 1)

    log_stream.info(' ')
    log_stream.info(' ==> ' + alg_name + ' (Version: ' + alg_version + ' Release_Date: ' + alg_release + ')')
    log_stream.info(' ==> TIME ELAPSED: ' + str(time_elapsed) + ' seconds')
    log_stream.info(' ==> ... END')
    log_stream.info(' ==> Bye, Bye')
    log_stream.info(' ============================================================================ ')
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to get script argument(s)
def get_args():
    parser_handle = ArgumentParser()
    parser_handle.add_argument('-settings_file', action="store", dest="alg_settings")
    parser_handle.add_argument('-time', action="store", dest="alg_time")
    parser_values = parser_handle.parse_args()

    alg_settings = 'configuration.json'
    if parser_values.alg_settings:
        alg_settings = parser_values.alg_settings

    alg_time = None
    if parser_values.alg_time:
        alg_time = parser_values.alg_time

    return alg_settings, alg_time

# -------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# Call script from external library
if __name__ == '__main__':
    main()
# ----------------------------------------------------------------------------
