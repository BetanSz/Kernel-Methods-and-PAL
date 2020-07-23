# --*-- coding:utf-8 --*--
"""
This module contain all FP-style functions.
FP: focuses on defining what to do, instead of how (that is left to the compiler which is more efficient than me).

Important function are here and it should be properly tested

"""
from tools_SG_common import *


def key_metric_data(data, metrics=None):
    """
    Returns a dictionary as expected by the class
    """
    return {name: metric(data) for name, metric in metrics.iteritems()}


def func_composition(*functions):
    def compose2(f, g):
        return lambda x: f(g(x))
    return reduce(compose2, functions, lambda x: x)


def cmpgtau(tau, f, c, w):
    """
    Descp
    """
    # this function can be modified or eliminated to work with  to_pppack_style
    Ni = [len(t) for t in tau]
    N = np.prod(Ni)
    gtau = np.zeros((N,), order='F')
    for x, i in zip(itertools.product(*tau), itertools.product(*[range(nx) for nx in Ni])):
        gtau[getidx(i, Ni)] = fgentz(x, f, c, w)
    return gtau


def getidx(i, Ni):
    """
    Get the position of the tuple *i* in the equivalent one-dimensional
    vector. The list *Ni* contains the maximum values of the indices,
    which start all from 0.
    """
    k, pos = len(Ni), 0
    for j in range(k - 1, 0, -1):
        pos = Ni[j - 1] * (i[j] + pos)
    return pos + i[0]


def normalize(vec, a_ext=None, b_ext=None):
    """
    normilizes vector
    unit-test
    a_ext<vec[0]<vec[-1]<b_ext
    reconvert to original vector and have the same values
    critical function
    """

    if a_ext == None and b_ext == None:
        a = vec[0]
        b = vec[-1]
    else:
        if not (a_ext < vec[0]) and (vec[0] <= vec[-1]) and (vec[-1] < b_ext):
            raise ValueError('extrema of real interval not valid')
        a = a_ext
        b = b_ext

    return [a, b], [(x - a) / (b - a) for x in vec]


def compute_error(f, g, *kwargs):
    """Computation of the error. Any difference of length will be filled with None's that will blow up in healthy way later"""
    # Order 1 actions
    # esto tiene q ser unificado forall unidad y tipo
    g = list(g)  # list comprehension, transofmration to list, tee for multiple generator compite for performance. Which is better requires testing
    # print type(f), type(g)
    if type(f) != type(g) or type(f) != list:
        raise ValueError('In the implementation today lists are required')

    if not kwargs:
        raise ValueError('no units given')
    # print kwargs
    # This could be done only for 'e' and then store a lazy eval
    ret_dict = OrderedDict()
    #ret_dict['results'] = OrderedDict()
    for key, unit in kwargs:
        ret_dict.update({key: [unit(fi, gi)
                               for fi, gi in itertools.izip_longest(f, g, fillvalue=None)]})
    return ret_dict


def cut_condition(f_i, g_i):
    """
    If this conditions are met the error calculations exists with an error
    """
    if abs(100 * (f_i - g_i) / g_i > 5) and 100 * f_i / g_i > 99 and 100 * f_i / g_i < 101:
        raise ValueError('analize excluded xs')
    if g_i <= 0:
        raise ValueError('Negative XS in statistics')
    return True

