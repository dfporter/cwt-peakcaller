"""
http://statsmodels.sourceforge.net/devel/examples/generated/example_gmle.html
"""
import numpy as np
from scipy.stats import nbinom
import scipy.stats
from statsmodels.base.model import GenericLikelihoodModel

import sys
import numpy.random
import scipy.stats
import time
import threading
import timeit
import scipy.optimize

def info(arr):
    return "mean: {mn} median: {med} variance: {var}".format(
        mn=np.nanmean(arr), med=np.nanmedian(arr), var=np.nanvar(arr)
    )

def likelihood_f(xxx_todo_changeme, x):
        (n,p) = xxx_todo_changeme
        return -1 * np.sum([scipy.stats.nbinom.logpmf(an_obs, n, p) for an_obs in x])

def fit_params(obs,just_positives=False, raise_background=False):
    if len(obs) == 0: return None
    obs = [int(x * 1e2) for x in obs]# if x > 0.]
    try:
        if np.max(obs) <= 0: return None
        pos_vals = [x for x in obs if (x > 0)]
    except:
        print("failure on {o}".format(o=obs))
        return None
    if just_positives:
        obs = [x for x in obs if x>0]
    elif raise_background:
        min_pos_obs = np.min(pos_vals)
        obs = [np.max([min_pos_obs, x]) for x in obs]
    rs = scipy.optimize.fmin(
        likelihood_f,
        x0=np.array((3, 0.5)),
        args=(obs,),maxfun=1000)
    return rs


def generate_data(n):
    out = []
    for i in range(n):
        out.append(
            scipy.stats.nbinom.rvs(1e3, 0.7, size=21)
        )
    return out

def process(_arr2d):
    params = []
    for r in _arr2d:
        params.append(fit_params(r))
    return params

def wrapper(process, *args, **kwargs):
    def wrapped():
        return process(*args, **kwargs)
    return wrapped


if __name__ == '__main__':
    num_arr=1
    arr2d = generate_data(10.)
    wrapped = wrapper(process, arr2d)
    took = timeit.timeit(wrapped, number=num_arr)/60.
    print("Time to fit {v} arrays: \
    {t} m. Would take {m} h to run on 6 * 1e4 arrays.".format(
        t=took, m=float(6*1e4*took)/(60.*num_arr*len(arr2d)),
        v=num_arr*len(arr2d)))
    


# .rvs gives a special (?) meaning to size, it is just the number of variates
# to generate.
# scale=1 by default. Called the scale parameter.
# loc=0 by default.
# vals = vals * scale + loc.
# Arguments passed to nbinom.rvs appear to go through rv_generic objects,
# which applies the loc and scale variables, regardless of the specific
# distribution.
