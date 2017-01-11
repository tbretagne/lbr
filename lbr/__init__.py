# This Python file uses the following encoding: utf-8

import lets_be_rational as lb
import normaldistribution as nd


def cum_norm(x):
    """ The cumulative normal distribution for the given argument x."""
    return nd.norm_cdf(x)


def inv_cum_corm(p):
    """ The inverse cumulative normal distribution for the given probability p. """
    return nd.inverse_norm_cdf(p)


def normalised_black_call(x, s):
    """ The normalised Black call option value [exp(x/2)·Phi(x/s+s/2)-exp(-x/2)·Phi(x/s-s/2)] with x=ln(F/K) and s=sigma·sqrt(T). """
    return lb.normalised_black_call(x, s)


def normalised_black(x, s, q):
    """ The normalised Black option value q·[exp(x/2)·Phi(q·(x/s+s/2))-exp(-x/2)·Phi(q·(x/s-s/2))] with x=ln(F/K) and s=sigma·sqrt(T). """
    """ q = ±1 for calls/puts"""
    return lb.normalised_black(x, s, q)


def black(F, K, sigma, T, q):
    """ The Black option value Black(F,K,sigma,T,q).
    Forward F, Strike K,
    Volatility sigma, Time Fraction T,
    q = ±1 """
    return lb.black(F, K, sigma, T, q)


def implied_volatility(price, F, K, T, q):
    """ The implied volatility s such that the given price equals the normalised Black option value [F·Phi(q·(x/s+s/2))-K·Phi(q·(x/s-s/2))] with x=ln(F/K) and s=sigma·sqrt(T). """
    return lb.implied_volatility_from_a_transformed_rational_guess(price, F, K, T, q)


def normalised_implied_volatility(beta, x, q):
    """ The implied volatility s such that the given value beta equals the normalised Black option value q·[exp(x/2)·Phi(q·(x/s+s/2))-exp(-x/2)·Phi(q·(x/s-s/2))] with x=ln(F/K) and s=sigma·sqrt(T). """
    """ Inputs:
        beta: Black/sqrt(F*K)
        x = ln(F/K)
        q = ±1 for calls/puts
    """
    return lb.normalised_implied_volatility_from_a_transformed_rational_guess(beta, x, q)
