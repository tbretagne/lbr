"""
 This source code resides at www.jaeckel.org/LetsBeRational.7z .

 ======================================================================================
 Copyright © 2013-2014 Peter Jäckel.

 Permission to use, copy, modify, and distribute this software is freely granted,
 provided that this notice is preserved.

 WARRANTY DISCLAIMER
 The Software is provided "as is" without warranty of any kind, either express or implied,
 including without limitation any implied warranties of condition, uninterrupted use,
 merchantability, fitness for a particular purpose, or non-infringement.
 ======================================================================================
"""
import numpy as np

DBL_EPSILON = 2.2204460492503131e-16
DBL_MAX = 1.79769e+308
DBL_MIN = -DBL_MAX

minimum_rational_cubic_control_parameter_value = -(1 - np.sqrt(DBL_EPSILON))
maximum_rational_cubic_control_parameter_value = 2 / (DBL_EPSILON * DBL_EPSILON)


def is_zero(x):

    return np.fabs(x) < DBL_MIN


def rational_cubic_interpolation(x, x_l, x_r, y_l, y_r, d_l, d_r, r):

    h = (x_r - x_l)
    if np.fabs(h) <= 0:
        return 0.5 * (y_l + y_r)
    # r should be greater than -1. We do not use  assert(r > -1)  here in order to allow values such as NaN to be propagated as they should.
    t = (x - x_l) / h
    if (not (r >= maximum_rational_cubic_control_parameter_value)):
        t = (x - x_l) / h
        omt = 1 - t
        t2 = t * t
        omt2 = omt * omt
    # Formula (2.4) divided by formula (2.5)
        return (y_r * t2 * t + (r * y_r - h * d_r) * t2 * omt + (r * y_l + h * d_l) * t * omt2 + y_l * omt2 * omt) / (1 + (r - 3) * t * omt)

    # Linear interpolation without over-or underflow.
    else:
        return y_r * t + y_l * (1 - t)


def rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_l):

    h = (x_r - x_l)
    numerator = 0.5 * h * second_derivative_l + (d_r - d_l)
    if is_zero(numerator):
        return 0
    else:
        denominator = (y_r - y_l) / h - d_l
        if is_zero(denominator):
            if numerator > 0:
                return maximum_rational_cubic_control_parameter_value
            else:
                return minimum_rational_cubic_control_parameter_value
        else:
            return numerator / denominator


def rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_r):

    h = (x_r - x_l)
    numerator = 0.5 * h * second_derivative_r + (d_r - d_l)
    if is_zero(numerator):
        return 0
    else:
        denominator = d_r - (y_r - y_l) / h
        if is_zero(denominator):
            if numerator > 0:
                return maximum_rational_cubic_control_parameter_value
            else:
                return minimum_rational_cubic_control_parameter_value
        else:
            return numerator / denominator


def minimum_rational_cubic_control_parameter(d_l, d_r, s, preferShapePreservationOverSmoothness):

    monotonic = (d_l * s >= 0) and (d_r * s >= 0)
    convex = (d_l <= s) and (s <= d_r)
    concave = (d_l >= s) and (s >= d_r)

    if (not monotonic and not convex and not concave):
        # If 3==r_non_shape_preserving_target, this means revert to standard cubic.
        return minimum_rational_cubic_control_parameter_value
    else:
        d_r_m_d_l = d_r - d_l
        d_r_m_s = d_r - s
        s_m_d_l = s - d_l
        r1 = -DBL_MAX
        r2 = r1
        # If monotonicity on this interval is possible, set r1 to satisfy the monotonicity condition (3.8).
        if monotonic:
            if not is_zero(s):  # 3.8), avoiding division by zero.
                r1 = (d_r + d_l) / s  # (3.8)
            elif preferShapePreservationOverSmoothness:
                # If division by zero would occur, and shape preservation is preferred, set value to enforce linear interpolation.
                r1 = maximum_rational_cubic_control_parameter_value  # This value enforces linear interpolation.

        if (convex or concave):
            if ((not is_zero(s_m_d_l)) or is_zero(d_r_m_s)):  # (3.18), avoiding division by zero.
                r2 = max(np.fabs(d_r_m_d_l / d_r_m_s),
                            np.fabs(d_r_m_d_l / s_m_d_l))
            elif preferShapePreservationOverSmoothness:
                r2 = maximum_rational_cubic_control_parameter_value  # This value enforces linear interpolation.
        elif (monotonic and preferShapePreservationOverSmoothness):
            r2 = maximum_rational_cubic_control_parameter_value  # This enforces linear interpolation along segments that are inconsistent with the slopes on the boundaries, e.g., a perfectly horizontal segment that has negative slopes on either edge.

        return max(minimum_rational_cubic_control_parameter_value,
                      max(r1, r2))


def convex_rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_l, preferShapePreservationOverSmoothness):

    r = rational_cubic_control_parameter_to_fit_second_derivative_at_left_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_l)
    r_min = minimum_rational_cubic_control_parameter(d_l, d_r, (y_r-y_l)/(x_r-x_l), preferShapePreservationOverSmoothness)

    return max(r, r_min)


def convex_rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_r, preferShapePreservationOverSmoothness):

    r = rational_cubic_control_parameter_to_fit_second_derivative_at_right_side(x_l, x_r, y_l, y_r, d_l, d_r, second_derivative_r)
    r_min = minimum_rational_cubic_control_parameter(d_l, d_r, (y_r-y_l)/(x_r-x_l), preferShapePreservationOverSmoothness)

    return max(r, r_min)
