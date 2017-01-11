
"""

 Original Fortran code taken from http:#www.netlib.org/specfun/erf, adapted by hand.
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
from erf_cody import erfc_cody
# from scipy.norm import pdf, ppf, cdf

DBL_EPSILON = 2.2204460492503131e-16
DBL_MAX = 1.79769e+308
DBL_MIN = -DBL_MAX

ONE_OVER_SQRT_TWO = 0.7071067811865475244008443621048490392848359376887
ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399460599343818684758586311649
SQRT_TWO_PI = 2.506628274631000502415765284811045253006986740610

norm_cdf_asymptotic_expansion_first_threshold = -10.0
norm_cdf_asymptotic_expansion_second_threshold = -1/np.sqrt(DBL_EPSILON)


def norm_pdf(x):

    return ONE_OVER_SQRT_TWO_PI * np.exp(-0.5 * x * x)


def norm_cdf(z):
    if z <= norm_cdf_asymptotic_expansion_first_threshold:
        sum_ = 1
        if z >= norm_cdf_asymptotic_expansion_second_threshold:
            zsqr = z * z
            i = 1
            g = 1
            a = DBL_MAX
            lasta = a
            while ((lasta > a) and (a >= np.fabs(sum_ * DBL_EPSILON))):
                lasta = a
                x = (4 * i - 3) / zsqr
                y = x * ((4 * i - 1) / zsqr)
                a = g * (x - y)
                sum_ -= a
                g *= y
                i += 1
                a = np.fabs(a)
        return -norm_pdf(z) * sum_ / z
    return 0.5*erfc_cody(-z*ONE_OVER_SQRT_TWO)


def inverse_norm_cdf(u):

    split1 = 0.425
    split2 = 5.0
    const1 = 0.180625
    const2 = 1.6

    # Coefficients for P close to 0.5
    A0 = 3.3871328727963666080E0
    A1 = 1.3314166789178437745E+2
    A2 = 1.9715909503065514427E+3
    A3 = 1.3731693765509461125E+4
    A4 = 4.5921953931549871457E+4
    A5 = 6.7265770927008700853E+4
    A6 = 3.3430575583588128105E+4
    A7 = 2.5090809287301226727E+3
    B1 = 4.2313330701600911252E+1
    B2 = 6.8718700749205790830E+2
    B3 = 5.3941960214247511077E+3
    B4 = 2.1213794301586595867E+4
    B5 = 3.9307895800092710610E+4
    B6 = 2.8729085735721942674E+4
    B7 = 5.2264952788528545610E+3
    # Coefficients for P not close to 0, 0.5 or 1.
    C0 = 1.42343711074968357734E0
    C1 = 4.63033784615654529590E0
    C2 = 5.76949722146069140550E0
    C3 = 3.64784832476320460504E0
    C4 = 1.27045825245236838258E0
    C5 = 2.41780725177450611770E-1
    C6 = 2.27238449892691845833E-2
    C7 = 7.74545014278341407640E-4
    D1 = 2.05319162663775882187E0
    D2 = 1.67638483018380384940E0
    D3 = 6.89767334985100004550E-1
    D4 = 1.48103976427480074590E-1
    D5 = 1.51986665636164571966E-2
    D6 = 5.47593808499534494600E-4
    D7 = 1.05075007164441684324E-9
    # Coefficients for P very close to 0 or 1
    E0 = 6.65790464350110377720E0
    E1 = 5.46378491116411436990E0
    E2 = 1.78482653991729133580E0
    E3 = 2.96560571828504891230E-1
    E4 = 2.65321895265761230930E-2
    E5 = 1.24266094738807843860E-3
    E6 = 2.71155556874348757815E-5
    E7 = 2.01033439929228813265E-7
    F1 = 5.99832206555887937690E-1
    F2 = 1.36929880922735805310E-1
    F3 = 1.48753612908506148525E-2
    F4 = 7.86869131145613259100E-4
    F5 = 1.84631831751005468180E-5
    F6 = 1.42151175831644588870E-7
    F7 = 2.04426310338993978564E-15

    if u <= 0:
        return np.log(u)
    elif u >= 1:
        return np.log(1 - u)
    else:
        q = u - 0.5
        if np.fabs(q) <= split1:
            r = const1 - q * q
            ret1 = q * (((((((A7 * r + A6) * r + A5) * r + A4) * r + A3) * r + A2) * r + A1) * r + A0)
            ret2 = (((((((B7 * r + B6) * r + B5) * r + B4) * r + B3) * r + B2) * r + B1) * r + 1.0)
            return ret1 / ret2
        else:
            if q < 0.:
                r = u
            else:
                r = 1. - u
            if r < split2:
                r = r - const2
                ret1 = (((((((C7 * r + C6) * r + C5) * r + C4) * r + C3) * r + C2) * r + C1) * r + C0)
                ret2 = (((((((D7 * r + D6) * r + D5) * r + D4) * r + D3) * r + D2) * r + D1) * r + 1.0)
                ret = ret1 / ret2
            else:
                r = r - split2
                ret1 = (((((((E7 * r + E6) * r + E5) * r + E4) * r + E3) * r + E2) * r + E1) * r + E0)
                ret2 = (((((((F7 * r + F6) * r + F5) * r + F4) * r + F3) * r + F2) * r + F1) * r + 1.0)
                ret = ret1 / ret2

            if q < 0:
                return -ret
            else:
                return ret
