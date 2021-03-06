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
