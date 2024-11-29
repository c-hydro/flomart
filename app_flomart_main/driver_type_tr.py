"""
Class Features

Name:          driver_type_tr
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20231114'
Version:       '1.0.0'
"""

# ----------------------------------------------------------------------------------------------------------------------
# libraries
import logging
import numpy as np

from lib_utils_tr import cmp_tr_general, cmp_tr_series
from lib_info_args import logger_name, time_format_algorithm

# logging
log_stream = logging.getLogger(logger_name)
# debugging
# import matplotlib.pylab as plt
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
# class driver type
class DriverType:
    # ------------------------------------------------------------------------------------------------------------------
    # empty class initialization
    section_description = None
    section_discharge_idx = None

    section_scenario_tr_corr_factor = None
    section_scenario_tr_par_a = None
    section_scenario_tr_par_b = None
    section_scenario_tr_par_tr = None
    section_scenario_tr_par_qr = None
    section_scenario_tr_min, section_scenario_tr_max = None, None

    section_scenario_q_t = None
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # initialization class
    def __init__(self,
                 method_scenario_tr='method_regional',
                 section_name='section_data',
                 section_scenario_tr_min=None, section_scenario_tr_max=None,
                 scenario_index_tag='scenario_index',
                 scenario_index_right_tag='scenario_index_right', scenario_index_left_tag='scenario_index_left',
                 scenario_weight_right_tag='scenario_weight_right', scenario_weight_left_tag='scenario_weight_left',
                 scenario_discharge_tag='scenario_discharge',
                 scenario_type_tag='scenario_type', scenario_time_tag='scenario_time',
                 scenario_n_tag='scenario_n', scenario_attrs_tag='scenario_attrs'):

        self.method_scenario_tr = method_scenario_tr
        self.section_name = section_name
        self.section_scenario_tr_min, self.section_scenario_tr_max = section_scenario_tr_min, section_scenario_tr_max

        self.scenario_index_tag = scenario_index_tag
        self.scenario_index_right_tag = scenario_index_right_tag
        self.scenario_index_left_tag = scenario_index_left_tag
        self.scenario_weight_right_tag = scenario_weight_right_tag
        self.scenario_weight_left_tag = scenario_weight_left_tag
        self.scenario_discharge_tag = scenario_discharge_tag
        self.scenario_type_tag = scenario_type_tag
        self.scenario_time_tag = scenario_time_tag
        self.scenario_n_tag = scenario_n_tag
        self.scenario_attrs_tag = scenario_attrs_tag

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # init scenario tr method
    def init_scenario_tr(self, section_data):

        # select method and variable(s)
        if self.method_scenario_tr == 'method_regional':
            # generic information
            self.section_description = section_data['section_description']
            self.section_discharge_idx = section_data['section_discharge_idx']
            # method information
            self.section_scenario_tr_corr_factor = section_data['section_par_correction_factor']
            self.section_scenario_tr_par_a = section_data['section_par_a']
            self.section_scenario_tr_par_b = section_data['section_par_a']
            self.section_scenario_tr_par_tr = section_data['section_par_tr']
            self.section_scenario_tr_par_qr = section_data['section_par_qr']
        elif self.method_scenario_tr == 'method_q_t':
            # generic information
            self.section_description = section_data['section_description']
            self.section_discharge_idx = section_data['section_discharge_idx']
            # method information
            self.section_scenario_q_t = section_data['section_q_t']
        else:
            log_stream.error(' ===> TR method "' + self.method_scenario_tr + '" is not allowed')
            raise NotImplemented('Case not implemented yet')

        # check section name
        assert self.section_description == self.section_name, ' ===> Section name is not the same, check the settings'
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to collect scenario tr method
    def collect_scenario_tr(self,
                            section_scenario_trs, section_scenario_trs_right, section_scenario_trs_left,
                            section_scenario_weight_right, section_scenario_weight_left,
                            section_discharge_values, section_type_values, section_discharge_times,
                            section_n_values, section_discharge_attrs):

        scenario_tr_obj = {
            self.scenario_index_tag: section_scenario_trs,
            self.scenario_index_right_tag: section_scenario_trs_right,
            self.scenario_index_left_tag: section_scenario_trs_left,
            self.scenario_weight_right_tag: section_scenario_weight_right,
            self.scenario_weight_left_tag: section_scenario_weight_left,
            self.scenario_discharge_tag: section_discharge_values,
            self.scenario_type_tag: section_type_values,
            self.scenario_time_tag: section_discharge_times,
            self.scenario_n_tag: section_n_values,
            self.scenario_attrs_tag: section_discharge_attrs}

        return scenario_tr_obj
    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to wrap scenario tr method (time and discharge)
    def wrap_scenario_tr(self, section_discharge_times, section_discharge_values):

        # info of discharge index and discharge values
        log_stream.info(' ------> QIndex: ' + str(self.section_discharge_idx) + ' [m^3/s]')
        log_stream.info(' ------> QMax: ' + str(np.max(section_discharge_values)) + ' [m^3/s]')

        # choose computing scenario tr method
        if self.method_scenario_tr == 'method_regional':

            # execute scenario tr based on fx
            section_scenario_trs, section_scenario_trs_right, section_scenario_trs_left, \
                section_scenario_weights_right, section_scenario_weights_left = \
                self.compute_scenario_tr_fx(
                    self.section_discharge_idx, section_discharge_times, section_discharge_values,
                    section_scenario_tr_corr_factor=self.section_scenario_tr_corr_factor,
                    section_scenario_tr_min=self.section_scenario_tr_min,
                    section_scenario_tr_max=self.section_scenario_tr_max)

        elif self.method_scenario_tr == 'method_q_t':

            # execute scenario tr based on series (q-t)
            section_scenario_trs, section_scenario_trs_right, section_scenario_trs_left, \
                section_scenario_weights_right, section_scenario_weights_left = \
                self.compute_scenario_tr_series(
                    section_discharge_times, section_discharge_values,
                    section_scenario_tr_min=self.section_scenario_tr_min,
                    section_scenario_tr_max=self.section_scenario_tr_max)

        else:
            log_stream.error(' ===> TR method "' + self.method_scenario_tr + '" is not allowed')
            raise NotImplemented('Case not implemented yet')

        # method to check scenario tr, tr_right and tr_left
        self.check_scenario_tr(section_scenario_trs, section_scenario_trs_right, section_scenario_trs_left)

        return (section_scenario_trs, section_scenario_trs_right, section_scenario_trs_left,
                section_scenario_weights_right, section_scenario_weights_left)

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # check scenario tr
    @staticmethod
    def check_scenario_tr(scenario_trs, scenario_trs_right, scenario_trs_left):
        # iterate over tr steps
        for id_tr, (trs, trs_right, trs_left ) in enumerate(zip(scenario_trs, scenario_trs_right, scenario_trs_left)):
            if trs_left > trs_right:
                logging.error(' ===> "tr_left" is greater than "tr_right" at step "' + str(id_tr) + '"')
                raise RuntimeError('The left side must be less than the right side')
            if (trs >= trs_left) and (trs <= trs_right):
                pass
            else:
                logging.error(' ===> "tr" must be greater than "tr_min" and less than "tr_max".'
                              ' At step "' + str(id_tr) + '" this condition is not satisfied')
                raise RuntimeError('Check the "tr", "tr_left" and "tr_right" values')

    # ------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to compute scenario tr
    def compute_scenario_tr_series(
            self, section_discharge_times, section_discharge_values,
            section_scenario_tr_min=1, section_scenario_tr_max=500):

        # check format of values and times
        if not isinstance(section_discharge_values, list):
            section_discharge_values = [section_discharge_values]
        if not isinstance(section_discharge_times, list):
            section_discharge_times = [section_discharge_times]

        # iterate over time and discharge values
        section_scenario_trs = []
        section_scenario_trs_right, section_scenario_trs_left = [], []
        section_scenario_weights_right, section_scenario_weights_left = [], []
        section_scenario_tr_check, section_scenario_tr_right_check, section_scenario_tr_left_check = [], [], []
        for section_discharge_id, (section_discharge_time, section_discharge_value) in enumerate(
                zip(section_discharge_times, section_discharge_values)):

            # info time and discharge start
            log_stream.info(' -------> Time: ' + str(section_discharge_time) +
                            ' :: Q(t): ' + str(section_discharge_value) + ' ... ')

            # check discharge value
            if section_discharge_value >= 0.0:

                # compute tr and weights
                (section_scenario_tr_rounded, section_scenario_tr_left, section_scenario_tr_right,
                 section_scenario_weight_left, section_scenario_weight_right) = cmp_tr_series(
                    section_discharge_value, self.section_scenario_q_t)

                # convert tr values to int (just in case)
                section_scenario_tr = int(section_scenario_tr_rounded)
                section_scenario_tr_right = int(section_scenario_tr_right)
                section_scenario_tr_left = int(section_scenario_tr_left)
                # check the tr values
                if section_scenario_tr_left > section_scenario_tr_right:
                    log_stream.error(' ===> "tr_left" is greater than "tr_right"')
                    raise RuntimeError('Check your tr values')

                # info scenario tr defined by the algorithm
                log_stream.info(' -------> tr: ' + str(section_scenario_tr) +
                                ', tr_left: ' + str(section_scenario_tr_left) + ' (weight: ' +
                                str(section_scenario_weight_left) + ')' +
                                ', tr_right: ' + str(section_scenario_tr_right) + ' (weight: ' +
                                str(section_scenario_weight_right) + ')')

                # check scenarios boundaries
                if section_scenario_tr < section_scenario_tr_min:

                    # scenarios tr less than minimum tr
                    section_scenario_tr = section_scenario_tr_min
                    section_scenario_tr_check.append(section_discharge_time.strftime(time_format_algorithm))

                elif section_scenario_tr > section_scenario_tr_max:

                    # scenarios tr greater than maximum tr
                    log_stream.warning(
                        ' ===> At time "' + section_discharge_time.strftime(time_format_algorithm) +
                        '" find the "tr ' + str(section_scenario_tr) +
                        '" greater then "tr_max ' + str(section_scenario_tr_max) + '"')
                    log_stream.warning(' ===> Assign "tr" == "tr_max')
                    section_scenario_tr = section_scenario_tr_max

                # check weight scenarios boundaries
                if section_scenario_tr_left < section_scenario_tr_min:

                    # scenarios weights less than minimum tr
                    section_scenario_tr_left_check.append(
                        section_discharge_time.strftime(time_format_algorithm))
                    # assign default conditions
                    section_scenario_tr_right = section_scenario_tr_min
                    section_scenario_tr_left = section_scenario_tr_min
                    section_scenario_weight_right, section_scenario_weight_left = 0.5, 0.5

                elif section_scenario_tr_right > section_scenario_tr_max:

                    # scenarios weights greater than maximum tr
                    log_stream.warning(
                        ' ===> At time "' + section_discharge_time.strftime(time_format_algorithm) +
                        '" find the "tr_right ' + str(section_scenario_tr_right) +
                        '" greater then "tr_max ' + str(section_scenario_tr_max) + '"')
                    log_stream.warning(' ===> Assign "tr_left" == "tr_max" and "tr_right" == "tr_max"')
                    # assign default conditions
                    section_scenario_tr_left = section_scenario_tr_max
                    section_scenario_tr_right = section_scenario_tr_max
                    section_scenario_weight_right, section_scenario_weight_left = 0.5, 0.5

            else:
                # discharge value less than zero
                section_scenario_tr = np.nan
                section_scenario_tr_right, section_scenario_tr_left = np.nan, np.nan
                section_scenario_weight_right, section_scenario_weight_left = np.nan, np.nan

            # store information step by step
            section_scenario_trs.append(section_scenario_tr)
            section_scenario_trs_right.append(section_scenario_tr_right)
            section_scenario_trs_left.append(section_scenario_tr_left)
            section_scenario_weights_right.append(section_scenario_weight_right)
            section_scenario_weights_left.append(section_scenario_weight_left)

            # info time and discharge end
            log_stream.info(' -------> Time: ' + str(section_discharge_time) +
                            ' :: Q(t): ' + str(section_discharge_value) + ' ... DONE')

        # check scenarios boundaries (generic)
        if section_scenario_tr_check:
            section_scenario_tr_str = ', '.join(section_scenario_tr_check)
            log_stream.warning(' ===> At times "' + section_scenario_tr_str +
                               '" found the "tr" less then "tr_min ' + str(section_scenario_tr_min) + '"')
            log_stream.warning(' ===> Set the "tr" equal to "tr_min"')
        # check scenarios boundaries (left)
        if section_scenario_tr_left_check:
            section_scenario_tr_str = ', '.join(section_scenario_tr_left_check)
            log_stream.warning(' ===> At times "' + section_scenario_tr_str +
                               '" found the "tr_left" less then "tr_min ' + str(section_scenario_tr_min) + '"')
            log_stream.warning(' ===> Set the "tr_left" equal to "tr_min"')

        return section_scenario_trs, section_scenario_trs_right, section_scenario_trs_left, \
            section_scenario_weights_right, section_scenario_weights_left

    # ----------------------------------------------------------------------------------------------------------------------

    # ------------------------------------------------------------------------------------------------------------------
    # method to compute scenario tr
    def compute_scenario_tr_fx(
            self, section_discharge_idx, section_discharge_times, section_discharge_values,
            section_scenario_tr_corr_factor=None,
            section_scenario_tr_min=1, section_scenario_tr_max=500):

        # check format of values and times
        if not isinstance(section_discharge_values, list):
            section_discharge_values = [section_discharge_values]
        if not isinstance(section_discharge_times, list):
            section_discharge_times = [section_discharge_times]

        # check correction factor
        if section_scenario_tr_corr_factor is None:
            log_stream.warning(' ===> The "section_scenario_tr_corr_factor" is defined by NoneType; '
                               'this could create some errors in the execution of the code')

        # conditions to activate fx computation
        active_fx = False
        if section_discharge_idx > 0.0:
            active_fx = True
        elif (section_discharge_idx <= 0.0) and (
                (self.section_scenario_tr_par_tr is not None) and
                (self.section_scenario_tr_par_qr is not None)):
            active_fx = True
        else:
            log_stream.warning(' ===> The scenario computation is not activated due to the datasets information. '
                               'Information are not enough to compute the scenario or the case is not implemented yet')

        # check the fx activation
        if active_fx:

            # iterate over time and discharge values
            section_scenario_trs = []
            section_scenario_trs_right, section_scenario_trs_left = [], []
            section_scenario_weights_right, section_scenario_weights_left = [], []
            section_scenario_tr_check, section_scenario_tr_right_check, section_scenario_tr_left_check = [], [], []
            for section_discharge_id, (section_discharge_time, section_discharge_value) in enumerate(
                    zip(section_discharge_times, section_discharge_values)):

                # info time and discharge start
                log_stream.info(' -------> Time: ' + str(section_discharge_time) +
                                ' :: Q(t): ' + str(section_discharge_value) + ' ... ')

                # check discharge value
                if section_discharge_value >= 0.0:

                    # MATTEO: replaced two functions cmp_tr_exp and cmp_tr_lin into unique function:
                    section_scenario_tr_rounded, \
                        section_scenario_tr_right, section_scenario_tr_left, \
                        section_scenario_weight_right, section_scenario_weight_left = cmp_tr_general(
                            section_discharge_idx,
                            section_discharge_value,
                            section_tr_par_a=self.section_scenario_tr_par_a,
                            section_tr_par_b=self.section_scenario_tr_par_b,
                            section_tr_par_tr=self.section_scenario_tr_par_tr,
                            section_tr_par_qr=self.section_scenario_tr_par_qr,
                            section_tr_approx=3,
                            correction_discharge_factor=self.section_scenario_tr_corr_factor)

                    # convert tr values to int (just in case)
                    section_scenario_tr = int(section_scenario_tr_rounded)
                    section_scenario_tr_right = int(section_scenario_tr_right)
                    section_scenario_tr_left = int(section_scenario_tr_left)

                    # info scenario tr defined by the algorithm
                    log_stream.info(' -------> tr: ' + str(section_scenario_tr) +
                                    ', tr_left: ' + str(section_scenario_tr_left) + ' (weight: ' +
                                    str(section_scenario_weight_right) + ')' +
                                    ', tr_right: ' + str(section_scenario_tr_right) + ' (weight: ' +
                                    str(section_scenario_weight_left) + ')')

                    # check scenarios boundaries
                    if section_scenario_tr < section_scenario_tr_min:
                        # scenarios tr less than minimum tr
                        section_scenario_tr = section_scenario_tr_min
                        section_scenario_tr_check.append(section_discharge_time.strftime(time_format_algorithm))
                    elif section_scenario_tr > section_scenario_tr_max:
                        # scenarios tr greater than maximum tr
                        log_stream.warning(' ===> At time "' + section_discharge_time.strftime(time_format_algorithm) +
                                           '" find the "tr ' + str(section_scenario_tr) +
                                           '" greater then "tr_max ' + str(section_scenario_tr_max) + '"')
                        log_stream.warning(' ===> Assign "tr" == "tr_max')
                        section_scenario_tr = section_scenario_tr_max

                    # check weight scenarios boundaries
                    if section_scenario_tr_left < section_scenario_tr_min:
                        # scenarios weights less than minimum tr
                        section_scenario_tr_left_check.append(section_discharge_time.strftime(time_format_algorithm))
                        # assign default conditions
                        section_scenario_tr_right = section_scenario_tr_min
                        section_scenario_tr_left = section_scenario_tr_min
                        section_scenario_weight_right, section_scenario_weight_left = 0.5, 0.5

                    elif section_scenario_tr_right > section_scenario_tr_max:
                        # scenarios weights greater than maximum tr
                        log_stream.warning(' ===> At time "' + section_discharge_time.strftime(time_format_algorithm) +
                                           '" find the "tr_right ' + str(section_scenario_tr_right) +
                                           '" greater then "tr_max ' + str(section_scenario_tr_max) + '"')
                        log_stream.warning(' ===> Assign "tr_left" == "tr_max" and "tr_right" == "tr_max"')
                        # assign default conditions
                        section_scenario_tr_left = section_scenario_tr_max
                        section_scenario_tr_right = section_scenario_tr_max
                        section_scenario_weight_right, section_scenario_weight_left = 0.5, 0.5

                else:
                    # discharge value less than zero
                    section_scenario_tr = np.nan
                    section_scenario_tr_right, section_scenario_tr_left = np.nan, np.nan
                    section_scenario_weight_right, section_scenario_weight_left = np.nan, np.nan

                # store information step by step
                section_scenario_trs.append(section_scenario_tr)
                section_scenario_trs_right.append(section_scenario_tr_right)
                section_scenario_trs_left.append(section_scenario_tr_left)
                section_scenario_weights_right.append(section_scenario_weight_right)
                section_scenario_weights_left.append(section_scenario_weight_left)

                # info time and discharge end
                log_stream.info(' -------> Time: ' + str(section_discharge_time) +
                                ' :: Q(t): ' + str(section_discharge_value) + ' ... DONE')

            # check scenarios boundaries (generic)
            if section_scenario_tr_check:
                section_scenario_tr_str = ', '.join(section_scenario_tr_check)
                log_stream.warning(' ===> At times "' + section_scenario_tr_str +
                                   '" found the "tr" less then "tr_min ' + str(section_scenario_tr_min) + '"')
                log_stream.warning(' ===> Set the "tr" equal to "tr_min"')
            # check scenarios boundaries (left)
            if section_scenario_tr_left_check:
                section_scenario_tr_str = ', '.join(section_scenario_tr_left_check)
                log_stream.warning(' ===> At times "' + section_scenario_tr_str +
                                   '" found the "tr_left" less then "tr_min ' + str(section_scenario_tr_min) + '"')
                log_stream.warning(' ===> Set the "tr_left" equal to "tr_min"')
        else:
            # discharge idx less than zero
            section_scenario_trs = [np.nan] * section_discharge_values.__len__()
            section_scenario_trs_right = [np.nan] * section_discharge_values.__len__()
            section_scenario_trs_left = [np.nan] * section_discharge_values.__len__()
            section_scenario_weights_right = [np.nan] * section_discharge_values.__len__()
            section_scenario_weights_left = [np.nan] * section_discharge_values.__len__()

        return section_scenario_trs, section_scenario_trs_right, section_scenario_trs_left, \
            section_scenario_weights_right, section_scenario_weights_left

# ----------------------------------------------------------------------------------------------------------------------
