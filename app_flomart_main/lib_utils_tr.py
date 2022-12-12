"""
Library Features:

Name:          lib_utils_tr
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20220208'
Version:       '1.0.0'
"""

#######################################################################################
# Libraries
import logging
import numpy as np
from numpy.linalg import lstsq
from copy import deepcopy

from lib_info_args import logger_name
from scipy.interpolate import interp1d # MATTEO: added this for linear interpolation

# Logging
log_stream = logging.getLogger(logger_name)
#######################################################################################



# -------------------------------------------------------------------------------------
# Method to compute tr:
def cmp_tr_general(section_discharge_idx, section_discharge_value, domain_name, section_name, section_tr_approx=3):

    # A different procedure is used here for the specific domains.
    if domain_name == 'Foglia':
        # REGIONALIZZAZ. MARCHE:
        # par_a = -0.22753
        # par_b = 0.67353
        # new:
        # par_a = 0.0306
        # par_b = 0.2586
        # section_scenario_tr_tmp = np.exp(
        #     (section_discharge_idx * par_a + section_discharge_value) / (section_discharge_idx * par_b))

        # A relation T= f(Q) for each section:
        TT = [0,2,5,10,20,50,100,150,200,500]

        if section_name == 'Foglia_Bronzo':
            QQ = [0,120.1326, 181.3846, 256.1861, 352.8854, 419.2873, 426.1650, 435.9561, 441.4092, 507.7208]
            # par_a = -94.8631
            # par_b = 79.0748
            # Qlimit =  121 # (approximated)
        elif section_name == 'Foglia_FlomartTorreCotogna':
            QQ = [0,124.41, 199.022, 277.173, 402.642, 489.577, 530.205, 554.212, 572.15, 654.73]
            # par_a = -65.1845
            # par_b = 98.5691
            # Qlimit =  133.5074
        elif section_name == 'Foglia_FlomartLaBadia':
            QQ = [0,148.9, 267.8, 359.6, 497.2, 658.0, 757.2, 784.1, 813.9, 1005.1]
            # par_a = -28.3702
            # par_b = 154.2956
            # Qlimit =  135.3198
        elif section_name == 'Foglia_FlomartCasellaMontecchio':
            QQ = [0, 156.4, 283.4, 387.6, 532.4, 710.0, 819.4, 841.4, 880.2, 1104.0]
            # par_a = -20.9963
            # par_b = 169.2869
            # Qlimit =  138.3371
        elif section_name == 'Foglia_Montecchio':
            QQ = [0, 176.4, 328.6, 471.9, 631.7, 834, 1010.6, 1107.2, 1124.2, 1361.1]
            # par_a = 9.2097
            # par_b = 218.3909
            # Qlimit =  142.1674
        elif section_name == 'Foglia_Foce':
            QQ = [0, 188.5, 361.5, 505.9, 682.7, 904.6, 1080.3, 1231.6, 1270.8, 1516.6]
            # par_a = 28.7061
            # par_b = 245.115
            # Qlimit = 143  # approximated from 141.1947.

        log_stream.info(' --------------> Apply linear interpolation to compute tr ...')
        section_scenario_tr_tmp = interp1d(QQ, TT)(section_discharge_value)


        # if section_discharge_value >= Qlimit:
        #     log_stream.info(' --------------> Q >= Q* = ' + str(Qlimit) + ' --> Apply Exp formula to compute tr ...')
        #     section_scenario_tr_tmp = np.exp(
        #         (section_discharge_value - par_a) / (par_b))
        # else:
        #     log_stream.info(' --------------> Q < Q* = ' + str(Qlimit) + ' --> Apply Lin formula to compute tr ...')
        #     section_discharge_idx_min = 0.0
        #     section_discharge_idx_max = Qlimit
        #     section_scenario_tr_min = 0.0
        #     section_scenario_tr_max = np.exp((Qlimit -par_a) / (par_b))
        #     points_coords = [(section_discharge_idx_min, section_scenario_tr_min),
        #                      (section_discharge_idx_max, section_scenario_tr_max)]
        #     points_x, points_y = zip(*points_coords)
        #     line_a = np.vstack([points_x, np.ones(len(points_x))]).T
        #     line_m, line_c = lstsq(line_a, points_y)[0]
        #     section_scenario_tr_tmp = line_m * section_discharge_value + line_c


        section_scenario_tr_rounded, section_scenario_tr_right, section_scenario_tr_left, \
        section_scenario_weight_right, section_scenario_weight_left = cmp_tr_weights(
            section_scenario_tr_tmp, decimal_tr=section_tr_approx)


    elif domain_name == 'Chienti':
        # par_a = 0.0306
        # par_b = 0.2586
        # section_scenario_tr_tmp = np.exp(
        #     (section_discharge_idx * par_a + section_discharge_value) / (section_discharge_idx * par_b))

        # Relation T= f(Q) for each section:
        TT = [0,2,5,10,20,50,100,150,200,500]

        if section_name == 'Chienti_FlomartPolverina':
            QQ = [0, 123.68, 207.59, 276.59, 347.46, 523.75, 617.18, 658, 682.52, 777.66]
            # par_a = - 9.5042
            # par_b = 126.4921
            # Qlimit =  97.1819
        elif section_name == 'Chienti_Diga_Borgiano':
            QQ = [0, 128.61, 232.61, 309.34, 385.68, 558.97, 696.56, 739.97, 761.77, 835.17]
            # par_a = -11.2098
            # par_b = 139.3251
            # Qlimit =  107.7826
        elif section_name == 'Chienti_Diga_LeGrazie':
            QQ = [0, 175.7, 279.1, 479.3, 586.9, 675.8, 823.5, 869.0, 924.8, 1083.3]
            # par_a = -61.3058
            # par_b = 163.7679
            # Qlimit =  174.8211
        elif section_name == 'Chienti_FlomartLeGrazie':
            QQ = [0, 185.83, 286.79, 395.11, 492.46, 611.08, 706.69, 798.11, 849.04, 984.28]
            # par_a = -60.8448
            # par_b = 145.8813
            # Qlimit =  161.9620
        elif section_name == 'Chienti_Chienti1':
            QQ = [0, 191.5, 297, 410.4, 515.8, 644.4, 728.2, 843.9, 889, 1018.6]
            # par_a = -62.7266
            # par_b = 152.3858
            # Qlimit =  168.3524
        elif section_name == 'Chienti_FlomartPrimaFiastra':
            QQ = [0, 196.7, 307.4, 423.6, 532.4, 667.3, 776.1, 867.9, 919.1, 1047.2]
            # par_a = -65.3607
            # par_b = 157.7583
            # Qlimit =  174.7104
        elif section_name == 'Chienti_FlomartDopoFiastra':
            QQ = [0, 236, 390.6, 521.5, 641.8, 819.4, 957.1, 1040, 1139.2, 1348.9]
            # par_a = -66.3924
            # par_b = 199.0409
            # Qlimit =  204.3570
        elif section_name == 'Chienti_Chienti2':
            QQ = [0, 243.9, 403.4, 539, 663.6, 849.4, 992.8, 1067, 1180.6, 1394.8]
            # par_a = -68.7793
            # par_b = 205.7328
            # Qlimit =  211.3824
        elif section_name == 'Chienti_FlomartTorrione':
            QQ = [0, 250.8, 417.4, 562.4, 712.4, 887.1, 1001.6, 1074.5, 1230, 1436.4]
            # par_a = -80.6381
            # par_b = 209.7554
            # Qlimit =  226.0294
        elif section_name == 'Chienti_FlomartMonteCosaro':
            QQ = [0, 264.3, 447.1, 602.6, 774, 919.4, 1057.7, 1134.8, 1341.9, 1502.8]
            # par_a = -91.1979
            # par_b = 221.3025
            # Qlimit =  244.5931
        elif section_name == 'Chienti_FoceChienti':
            QQ = [0, 317, 549.4, 744.9, 931.2, 1163.4, 1363.6, 1425.1, 1461.4, 1643.2]
            # par_a = -172.8443
            # par_b = 246.4216
            # Qlimit = 343.6508

        log_stream.info(' --------------> Apply linear interpolation to compute tr ...')
        section_scenario_tr_tmp = interp1d(QQ, TT)(section_discharge_value)

        # if section_discharge_value >= Qlimit:
        #     log_stream.info(' --------------> Q >= Q* = ' + str(Qlimit) + ' --> Apply Exp formula to compute tr ...')
        #     section_scenario_tr_tmp = np.exp((par_a + section_discharge_value) / (par_b))
        # else:
        #     log_stream.info(' --------------> Q < Q* = ' + str(Qlimit) + ' --> Apply Lin formula to compute tr ...')
        #     section_discharge_idx_min = 0.0
        #     section_discharge_idx_max = Qlimit
        #     section_scenario_tr_min = 0.0
        #     section_scenario_tr_max = np.exp((par_a + Qlimit) / (par_b))
        #     points_coords = [(section_discharge_idx_min, section_scenario_tr_min),
        #                      (section_discharge_idx_max, section_scenario_tr_max)]
        #     points_x, points_y = zip(*points_coords)
        #     line_a = np.vstack([points_x, np.ones(len(points_x))]).T
        #     line_m, line_c = lstsq(line_a, points_y)[0]
        #     section_scenario_tr_tmp = line_m * section_discharge_value + line_c

        section_scenario_tr_rounded, section_scenario_tr_right, section_scenario_tr_left, \
        section_scenario_weight_right, section_scenario_weight_left = cmp_tr_weights(
            section_scenario_tr_tmp, decimal_tr=section_tr_approx)

    else:
        # for example:  domain_name == 'Entella'  REGIONALIZZAZ. LIGURIA:
        par_a = 0.5239
        par_b = 1.0433
        if section_discharge_value >= section_discharge_idx * self.correction_discharge_factor:
            log_stream.info(' --------------> Q >= Qindex * K = ' + str(
                section_discharge_idx * self.correction_discharge_factor) + ' --> Apply Exp formula to compute tr ...')
            section_scenario_tr_tmp = np.exp(
                (section_discharge_idx * par_a + section_discharge_value) / (section_discharge_idx * par_b))
        else:
            log_stream.info(' --------------> Q < Qindex * K = ' + str(
                section_discharge_idx * self.correction_discharge_factor) + ' --> Apply Lin formula to compute tr ...')
            section_scenario_tr_max = np.exp(
                (section_discharge_idx * par_a + section_discharge_idx * self.correction_discharge_factor) / (
                            section_discharge_idx * par_b))

            section_discharge_idx_min = 0.0
            section_discharge_idx_max = section_discharge_idx * self.correction_discharge_factor
            section_scenario_tr_min = 0.0
            points_coords = [(section_discharge_idx_min, section_scenario_tr_min),
                             (section_discharge_idx_max, section_scenario_tr_max)]
            points_x, points_y = zip(*points_coords)
            line_a = np.vstack([points_x, np.ones(len(points_x))]).T
            line_m, line_c = lstsq(line_a, points_y)[0]
            section_scenario_tr_tmp = line_m * section_discharge_value + line_c

        section_scenario_tr_rounded, section_scenario_tr_right, section_scenario_tr_left, \
        section_scenario_weight_right, section_scenario_weight_left = cmp_tr_weights(
            section_scenario_tr_tmp, decimal_tr=section_tr_approx)

    return section_scenario_tr_rounded, section_scenario_tr_right, section_scenario_tr_left,\
           section_scenario_weight_right, section_scenario_weight_left





# -------------------------------------------------------------------------------------
# Method to compute tr values and weights
def cmp_tr_weights(value_tr, decimal_tr=3):

    section_scenario_tr_cmp = deepcopy(value_tr)
    section_scenario_tr_left = np.ceil(section_scenario_tr_cmp)
    section_scenario_tr_right = np.floor(section_scenario_tr_cmp)

    section_scenario_weight_right = np.around(1 - (section_scenario_tr_cmp - section_scenario_tr_right), decimals=decimal_tr)
    section_scenario_weight_left = np.around(1 - (section_scenario_tr_left - section_scenario_tr_cmp), decimals=decimal_tr)

    section_scenario_tr_rounded = np.round(section_scenario_tr_cmp)

    return section_scenario_tr_rounded, \
        section_scenario_tr_right, section_scenario_tr_left, section_scenario_weight_right, section_scenario_weight_left
# -------------------------------------------------------------------------------------







# # -------------------------------------------------------------------------------------
# # Method to compute tr exp mode
# def cmp_tr_exp(section_discharge_idx, section_discharge_value, domain_name, section_name, section_tr_approx=3):
#
#     # section_scenario_tr = np.round(np.exp(
#     #    (section_discharge_idx * 0.5239 + section_discharge_value) / (section_discharge_idx * 1.0433)))
#
#
#
#     return section_scenario_tr_rounded, section_scenario_tr_right, section_scenario_tr_left, \
#         section_scenario_weight_right, section_scenario_weight_left
# # -------------------------------------------------------------------------------------





# -------------------------------------------------------------------------------------
# Method to compute tr linear mode
# def cmp_tr_linear(section_discharge_idx, section_discharge_value, domain_name,
#                   section_discharge_factor=1.16, section_tr_approx=3):
#
#
#     if domain_name[0] == 'Foglia':
#         # REGIONALIZZAZ. MARCHE:
#         # par_a = -0.22753
#         # par_b = 0.67353
#
#         # new:
#         # par_a = 0.0306
#         # par_b = 0.2586
#
#         if section_name == 'Foglia_Bronzo':
#             par_a = -94.8631
#             par_b = 79.0748
#         elif section_name == 'Foglia_FlomartTorreCotogna':
#             par_a = -65.1845
#             par_b = 98.5691
#         elif section_name == 'Foglia_FlomartLaBadia':
#             par_a = -28.3702
#             par_b = 154.2956
#         elif section_name == 'Foglia_FlomartCasellaMontecchio':
#             par_a = -20.9963
#             par_b = 169.2869
#         elif section_name == 'Foglia_Montecchio':
#             par_a = 9.2097
#             par_b = 218.3909
#         elif section_name == 'Foglia_Foce':
#             par_a = 28.7061
#             par_b = 245.115
#
#         section_scenario_tr_tmp = np.exp((par_a + section_discharge_value) / (par_b))
#         section_discharge_idx_min = 0.0
#         section_discharge_idx_max = section_discharge_factor
#
#     elif domain_name[0] == 'Chienti':
#         if section_name == 'Chienti_FlomartPolverina':
#             par_a = - 9.5042
#             par_b = 126.4921
#         elif section_name == 'Chienti_Diga_Borgiano':
#             par_a = -11.2098
#             par_b = 139.3251
#         elif section_name == 'Chienti_Diga_LeGrazie':
#             par_a = -61.3058
#             par_b = 163.7679
#         elif section_name == 'Chienti_FlomartLeGrazie':
#             par_a = -60.8448
#             par_b = 145.8813
#         elif section_name == 'Chienti_Chienti1':
#             par_a = -62.7266
#             par_b = 152.3858
#         elif section_name == 'Chienti_FlomartPrimaFiastra':
#             par_a = -65.3607
#             par_b = 157.7583
#         elif section_name == 'Chienti_FlomartDopoFiastra':
#             par_a = -66.3924
#             par_b = 199.0409
#         elif section_name == 'Chienti_Chienti2':
#             par_a = -68.7793
#             par_b = 205.7328
#         elif section_name == 'Chienti_FlomartTorrione':
#             par_a = -80.6381
#             par_b = 209.7554
#         elif section_name == 'Chienti_FlomartMonteCosaro':
#             par_a = -91.1979
#             par_b = 221.3025
#         elif section_name == 'Chienti_FoceChienti':
#             par_a = -172.8443
#             par_b = 246.4216
#
#         section_scenario_tr_tmp = np.exp((par_a + section_discharge_value) / (par_b))
#         section_discharge_idx_min = 0.0
#         section_discharge_idx_max = section_discharge_factor
#
#
#     elif domain_name[0] == 'Entella':
#         # REGIONALIZZAZ. LIGURIA:
#         # par_a = 0.5239
#         # par_b = 1.0433
#         section_scenario_tr_max = np.exp(
#             (section_discharge_idx * par_a + section_discharge_idx * section_discharge_factor) / (section_discharge_idx * par_b))
#
#         section_discharge_idx_min = 0.0
#         section_discharge_idx_max = section_discharge_idx * section_discharge_factor
#
#
#
#
#     section_scenario_tr_min = 0.0
#     points_coords = [(section_discharge_idx_min, section_scenario_tr_min),
#                      (section_discharge_idx_max, section_scenario_tr_max)]
#
#     points_x, points_y = zip(*points_coords)
#     line_a = np.vstack([points_x, np.ones(len(points_x))]).T
#
#     line_m, line_c = lstsq(line_a, points_y)[0]
#
#     section_scenario_tr_tmp = line_m * section_discharge_value + line_c
#
#     section_scenario_tr_rounded, section_scenario_tr_right, section_scenario_tr_left, \
#         section_scenario_weight_right, section_scenario_weight_left = cmp_tr_weights(
#             section_scenario_tr_tmp, decimal_tr=section_tr_approx)
#
#     return section_scenario_tr_rounded, \
#         section_scenario_tr_right, section_scenario_tr_left, \
#         section_scenario_weight_right, section_scenario_weight_left
# -------------------------------------------------------------------------------------



