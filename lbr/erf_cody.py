import numpy as np


def d_int(x):
    if x > 0:
        return np.floor(x)
    else:
        return -np.floor(-x)


def calerf(x, jint):

    a = [3.1611237438705656, 113.864154151050156, 377.485237685302021,
         3209.37758913846947, 0.185777706184603153]
    b = [23.6012909523441209, 244.024637934444173, 1282.61652607737228,
         2844.23683343917062]
    c__ = [0.564188496988670089, 8.88314979438837594, 66.1191906371416295,
           298.635138197400131, 881.95222124176909, 1712.04761263407058,
           2051.07837782607147, 1230.33935479799725, 2.15311535474403846e-8]
    d__ = [15.7449261107098347, 117.693950891312499, 537.181101862009858,
           1621.38957456669019, 3290.79923573345963, 4362.61909014324716,
           3439.36767414372164, 1230.33935480374942]
    p = [0.305326634961232344, 0.360344899949804439, 0.125781726111229246,
         0.0160837851487422766, 6.58749161529837803e-4, 0.0163153871373020978]
    q = [2.56852019228982242, 1.87295284992346047, 0.527905102951428412,
         0.0605183413124413191, 0.00233520497626869185]

    zero = 0.
    half = 0.5
    one = 1.
    two = 2.
    four = 4.
    sqrpi = 0.56418958354775628695
    thresh = 0.46875
    sixten = 16.

    xinf = 1.79e308
    xneg = -26.628
    xsmall = 1.11e-16
    xbig = 26.543
    xhuge = 6.71e7
    xmax = 2.53e307

    flag = 0

    y = np.fabs(x)

    if (y <= thresh):
        ysq = zero

        if (y > xsmall):
            ysq = y * y

        xnum = a[4] * ysq

        xden = ysq

        for i__ in range(1, 4):
            xnum = (xnum + a[i__ - 1]) * ysq
            xden = (xden + b[i__ - 1]) * ysq

        result = x * (xnum + a[3]) / (xden + b[3])

        if (jint != 0):
            result = one - result

        if (jint == 2):
            result = np.exp(ysq) * result

        flag = 800

    elif (y <= four):

        xnum = c__[8] * y

        xden = y

        for i__ in range(1, 8):

            xnum = (xnum + c__[i__ - 1]) * y

            xden = (xden + d__[i__ - 1]) * y

        result = (xnum + c__[7]) / (xden + d__[7])

        if (jint != 2):

            d__1 = y * sixten
            ysq = d_int(d__1) / sixten
            del__1 = (y - ysq) * (y + ysq)
            d__1 = np.exp(-ysq * ysq) * np.exp(-del__1)
            result = d__1 * result

    else:

        result = zero

        if (y >= xbig):

            if ((jint != 2) or (y >= xmax)):
                flag = 300

            if (y >= xhuge):

                result = sqrpi / y
                flag = 300

        if flag == 0:

            ysq = one / (y * y)

            xnum = p[5] * ysq
            xden = ysq

            for i__ in range(1, 5):

                xnum = (xnum + p[i__ - 1]) * ysq

                xden = (xden + q[i__ - 1]) * ysq

            result = ysq * (xnum + p[4]) / (xden + q[4])

            result = (sqrpi - result) / y

            if (jint != 2):

                d__1 = y * sixten
                ysq = d_int(d__1) / sixten

                del__1 = (y - ysq) * (y + ysq)

                d__1 = np.exp(-ysq * ysq) * np.exp(-del__1)
                result = d__1 * result

    if flag <= 300:
        if (jint == 0):
            result = (half - result) + half

        if (x < zero):
            result = -(result)

        elif (jint == 1):

            if (x < zero):
                result = two - result

        else:

            if (x < zero):

                if (x < xneg):

                    result = xinf
                else:
                    d__1 = x * sixten
                    ysq = d_int(d__1) / sixten

                    del__1 = (x - ysq) * (x + ysq)

                    y = np.exp(ysq * ysq) * np.exp(del__1)

                    result = y + y - result

    return result


def erf_cody(x):
    return calerf(x, 0)


def erfc_cody(x):
    return calerf(x, 1)


def erfcx_cody(x):
    return calerf(x, 2)
