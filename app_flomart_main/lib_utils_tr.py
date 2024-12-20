"""
Library Features:

Name:          lib_utils_tr
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231116'
Version:       '1.5.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import pandas as pd
import numpy as np

from numpy.linalg import lstsq
from copy import deepcopy

from lib_info_args import logger_name
from scipy.interpolate import interp1d # MATTEO: added this for linear interpolation

# Logging
log_stream = logging.getLogger(logger_name)
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to get tr parameters
def get_tr_params(section_name):

    # initialize parameters
    par_a, par_b, par_correction_factor, tr, qr = None, None, None, None, None

    # A relation T= f(Q) for each section:
    if section_name == 'Foglia_Bronzo':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 120.1326, 181.3846, 256.1861, 352.8854, 419.2873, 426.1650, 435.9561, 441.4092, 507.7208]
    elif section_name == 'Foglia_FlomartTorreCotogna':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 124.41, 199.022, 277.173, 402.642, 489.577, 530.205, 554.212, 572.15, 654.73]
    elif section_name == 'Foglia_FlomartLaBadia':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 148.9, 267.8, 359.6, 497.2, 658.0, 757.2, 784.1, 813.9, 1005.1]
    elif section_name == 'Foglia_FlomartCasellaMontecchio':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 156.4, 283.4, 387.6, 532.4, 710.0, 819.4, 841.4, 880.2, 1104.0]
    elif section_name == 'Foglia_Montecchio':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 176.4, 328.6, 471.9, 631.7, 834, 1010.6, 1107.2, 1124.2, 1361.1]
    elif section_name == 'Foglia_Foce':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 188.5, 361.5, 505.9, 682.7, 904.6, 1080.3, 1231.6, 1270.8, 1516.6]
    elif section_name == 'Chienti_FlomartPolverina':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 123.68, 207.59, 276.59, 347.46, 523.75, 617.18, 658, 682.52, 777.66]
    elif section_name == 'Chienti_Diga_Borgiano':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 128.61, 232.61, 309.34, 385.68, 558.97, 696.56, 739.97, 761.77, 835.17]
    elif section_name == 'Chienti_Diga_LeGrazie':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 175.7, 279.1, 479.3, 586.9, 675.8, 823.5, 869.0, 924.8, 1083.3]
    elif section_name == 'Chienti_FlomartLeGrazie':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 185.83, 286.79, 395.11, 492.46, 611.08, 706.69, 798.11, 849.04, 984.28]
    elif section_name == 'Chienti_Chienti1':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 191.5, 297, 410.4, 515.8, 644.4, 728.2, 843.9, 889, 1018.6]
    elif section_name == 'Chienti_FlomartPrimaFiastra':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 196.7, 307.4, 423.6, 532.4, 667.3, 776.1, 867.9, 919.1, 1047.2]
    elif section_name == 'Chienti_FlomartDopoFiastra':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 236, 390.6, 521.5, 641.8, 819.4, 957.1, 1040, 1139.2, 1348.9]
    elif section_name == 'Chienti_Chienti2':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 243.9, 403.4, 539, 663.6, 849.4, 992.8, 1067, 1180.6, 1394.8]
    elif section_name == 'Chienti_FlomartTorrione':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 250.8, 417.4, 562.4, 712.4, 887.1, 1001.6, 1074.5, 1230, 1436.4]
    elif section_name == 'Chienti_FlomartMonteCosaro':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 264.3, 447.1, 602.6, 774, 919.4, 1057.7, 1134.8, 1341.9, 1502.8]
    elif section_name == 'Chienti_FoceChienti':
        tr = [0, 2, 5, 10, 20, 50, 100, 150, 200, 500]
        qr = [0, 317, 549.4, 744.9, 931.2, 1163.4, 1363.6, 1425.1, 1461.4, 1643.2]
    elif section_name == 'Misa_Nevola':
        tr = [0, 5, 10, 20, 50, 100, 200, 500]
        qr = [0, 149, 200, 244,	313, 355, 393, 476]
    elif section_name == 'Misa_Misa':
        tr = [0, 5, 10, 20, 50, 100, 200, 500]
        qr = [0, 322, 424, 531, 663, 772, 882, 971]
    elif section_name == 'Misa_PonteGaribaldi':
        tr = [0, 5, 10, 20, 50, 100, 200, 500]
        qr = [0, 347, 449, 569,	722, 826, 949, 1052]
    elif section_name == 'Misa_FlomartCasine':
        tr = [0, 5, 10, 20, 50, 100, 200, 500]
        qr = [0, 311, 395, 456, 539, 576, 628, 735]
    elif section_name == 'Misa_FlomartBrugnetto':
        tr = [0, 5, 10, 20, 50, 100, 200, 500]
        qr = [0, 161, 210, 259, 327, 367, 406, 492]
    elif section_name == 'Misa_FlomartContradaMolino':
        tr = [0, 5, 10, 20, 50, 100, 200, 500]
        qr = [0, 289, 354, 416, 486, 526, 563, 656]
    else:
        par_a, par_b = 0.5239, 1.0433
        par_correction_factor = 1.16

    return par_a, par_b, par_correction_factor, tr, qr
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to compute tr values and weights based on series
def cmp_tr_series(section_discharge_value, section_q_t,
                  var_q='Q', var_t='T',
                  decimal_discharge=2, decimal_weight=3, decimal_tr=3):

    # organize discharge information
    section_discharge_min, section_discharge_max = section_q_t[var_q].min(), section_q_t[var_q].max()
    section_discharge_min = np.round(section_discharge_min, decimals=decimal_discharge)
    section_discharge_max = np.round(section_discharge_max, decimals=decimal_discharge)
    section_discharge_value = np.round(section_discharge_value, decimals=decimal_discharge)

    # compute tr and weight(s)
    if (section_discharge_value > section_discharge_min) and (section_discharge_value < section_discharge_max):

        section_qt_bounds = section_q_t.iloc[(section_q_t[var_q] - section_discharge_value).abs().argsort()[:2]]
        section_qt_bounds = section_qt_bounds.sort_values(var_q)

        section_discharge_lower = section_qt_bounds.iloc[0][var_q]
        section_scenario_tr_lower = section_qt_bounds.iloc[0][var_t]
        section_discharge_upper = section_qt_bounds.iloc[1][var_q]
        section_scenario_tr_upper = section_qt_bounds.iloc[1][var_t]

        section_scenario_weight_lower = np.abs(1 - np.abs(
            (section_discharge_value - section_discharge_lower) / (section_discharge_upper - section_discharge_lower)))
        section_scenario_weight_upper = np.abs(1 - np.abs(
            (section_discharge_value - section_discharge_upper) / (section_discharge_upper - section_discharge_lower)))

        section_scenario_tr_tmp = (section_scenario_tr_lower * section_scenario_weight_lower +
                                   section_scenario_tr_upper * section_scenario_weight_upper)

        section_scenario_tr_rounded = np.round(section_scenario_tr_tmp)

    elif section_discharge_value <= section_discharge_min:

        section_scenario_tr_min = section_q_t['T'].min()
        section_scenario_tr_rounded = section_scenario_tr_min
        section_scenario_tr_lower, section_scenario_tr_upper = section_scenario_tr_min, section_scenario_tr_min
        section_scenario_weight_lower, section_scenario_weight_upper = 0.5, 0.5

    elif section_discharge_value >= section_discharge_max:

        section_scenario_tr_max = section_q_t['T'].max()
        section_scenario_tr_rounded = section_scenario_tr_max
        section_scenario_tr_lower, section_scenario_tr_upper = section_scenario_tr_max, section_scenario_tr_max
        section_scenario_weight_lower, section_scenario_weight_upper = 0.5, 0.5

    else:
        log_stream.error(' ===> Discharge value is not supported by the procedure')
        raise NotImplemented('Case not implemented yet')

    # check the tr values boundaries
    if (section_scenario_tr_rounded >= section_scenario_tr_lower) and (section_scenario_tr_rounded <= section_scenario_tr_upper):
        pass
    else:
        logging.error(' ===> The tr value is not in the range of tr left and right')
        raise NotImplemented('Case not implemented yet')

    section_scenario_weight_lower = np.round(section_scenario_weight_lower, decimals=decimal_weight)
    section_scenario_weight_upper = np.round(section_scenario_weight_upper, decimals=decimal_weight)

    section_scenario_tr_rounded = np.round(section_scenario_tr_rounded, decimals=decimal_tr)
    section_scenario_tr_lower = np.round(section_scenario_tr_lower, decimals=decimal_tr)
    section_scenario_tr_upper = np.round(section_scenario_tr_upper, decimals=decimal_tr)

    return section_scenario_tr_rounded, \
        section_scenario_tr_lower, section_scenario_tr_upper, \
        section_scenario_weight_lower, section_scenario_weight_upper
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
# Method to compute tr:
def cmp_tr_general(section_discharge_idx, section_discharge_value,
                   section_tr_par_a=None, section_tr_par_b=None,
                   section_tr_par_tr=None, section_tr_par_qr=None,
                   section_tr_approx=3, correction_discharge_factor=None):

    if (section_tr_par_tr and section_tr_par_qr) and (section_tr_par_a is None and section_tr_par_b is None):

        # for example:  domain_name == 'Chienti, Foglia, Misa'  REGIONALIZZAZ. MARCHE

        log_stream.info(' --------> Apply linear interpolation to compute tr ... ')

        fx_scenario_tr = interp1d(section_tr_par_qr, section_tr_par_tr)
        section_scenario_tr_tmp = fx_scenario_tr(section_discharge_value)

        log_stream.info(' --------> Apply linear interpolation to compute tr ... DONE')

        section_scenario_tr_rounded, section_scenario_tr_right, section_scenario_tr_left, \
            section_scenario_weight_right, section_scenario_weight_left = cmp_tr_weights(
                section_scenario_tr_tmp, decimal_tr=section_tr_approx)

    elif (section_tr_par_tr is None and section_tr_par_qr is None) and (section_tr_par_a and section_tr_par_b):

        # for example:  domain_name == 'Entella'  REGIONALIZZAZ. LIGURIA
        # section_tr_par_a, section_tr_par_b = 0.5239, 1.0433

        # check of correction factor
        if correction_discharge_factor is None:
            log_stream.error(' ===> The "correction_discharge_factor" K is NoneType.')
            raise RuntimeError('The "correction_discharge_factor" must be defined by a value.')
        else:
            log_stream.info(' --------> The "correction_discharge_factor" K: ' + str(correction_discharge_factor))

        # case definition(s)
        if section_discharge_value >= section_discharge_idx * correction_discharge_factor:

            # apply linear relationship
            log_stream.info(' --------> Q >= Qindex * K = ' + str(
                section_discharge_idx * correction_discharge_factor) + ' :: Apply Exp formula to compute tr ... ')

            section_scenario_tr_tmp = np.exp(
                (section_discharge_idx * section_tr_par_a + section_discharge_value) /
                    (section_discharge_idx * section_tr_par_b))

            log_stream.info(' --------> Q >= Qindex * K = ' + str(
                section_discharge_idx * correction_discharge_factor) + ' :: Apply Exp formula to compute tr ... DONE')

        else:
            # apply exponential relationship
            log_stream.info(' --------> Q < Qindex * K = ' + str(
                section_discharge_idx * correction_discharge_factor) + ' :: Apply Lin formula to compute tr ... ')

            section_scenario_tr_max = np.exp(
                (section_discharge_idx * section_tr_par_a + section_discharge_idx * correction_discharge_factor) / (
                            section_discharge_idx * section_tr_par_b))

            section_discharge_idx_min = 0.0
            section_discharge_idx_max = section_discharge_idx * correction_discharge_factor
            section_scenario_tr_min = 0.0
            points_coords = [(section_discharge_idx_min, section_scenario_tr_min),
                             (section_discharge_idx_max, section_scenario_tr_max)]
            points_x, points_y = zip(*points_coords)
            line_a = np.vstack([points_x, np.ones(len(points_x))]).T
            line_m, line_c = lstsq(line_a, points_y)[0]
            section_scenario_tr_tmp = line_m * section_discharge_value + line_c

            log_stream.info(' --------> Q < Qindex * K = ' + str(
                section_discharge_idx * correction_discharge_factor) + ' :: Apply Lin formula to compute tr ... DONE')

        section_scenario_tr_rounded, section_scenario_tr_right, section_scenario_tr_left, \
            section_scenario_weight_right, section_scenario_weight_left = cmp_tr_weights(section_scenario_tr_tmp,
                                                                                         decimal_tr=section_tr_approx)

    else:
        log_stream.error(' ===> The arguments to compute the tr are not defined uniquely')
        raise RuntimeError('The parameters must be defined to select the method.')

    return section_scenario_tr_rounded, \
        section_scenario_tr_right, section_scenario_tr_left, \
        section_scenario_weight_right, section_scenario_weight_left

# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# Method to compute tr values and weights
def cmp_tr_weights(value_tr, decimal_tr=3):

    # define tr values (cmp, left, right)
    section_scenario_tr_cmp = deepcopy(value_tr)
    section_scenario_tr_left = np.floor(section_scenario_tr_cmp)    # lower side
    section_scenario_tr_right = np.ceil(section_scenario_tr_cmp)    # upper side

    # lower side
    section_scenario_weight_left = np.around(1 - (section_scenario_tr_cmp - section_scenario_tr_left),
                                             decimals=decimal_tr)
    # upper side
    section_scenario_weight_right = np.around(1 - (section_scenario_tr_right - section_scenario_tr_cmp),
                                              decimals=decimal_tr)

    # compute tr rounded
    section_scenario_tr_rounded = np.round(section_scenario_tr_cmp)

    return section_scenario_tr_rounded, \
        section_scenario_tr_right, section_scenario_tr_left, section_scenario_weight_right, section_scenario_weight_left
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# method to get tr values from ascii file
def get_tr_file(file_name, file_sep=';'):
    file_data = pd.read_csv(file_name, sep=file_sep)
    return file_data
# ----------------------------------------------------------------------------------------------------------------------
