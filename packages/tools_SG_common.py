# !/usr/bin/env python
"""
File containes common function to the entire script, toghether with common importation of modules.
Also many auxiliary functions.

Many function here are used in scripts that are not public including:
* Error analysis
* MPI application
* Plotting

"""

import os
import sys  # this package is un-recommended in literature
import time
import copy as cp
import functools
import multiprocessing as mp
import warnings
import traceback
import matplotlib.pyplot as plt
import matplotlib
import random
import pandas as pd
import dis
import toolz
import cytoolz as cz
from pprint import pprint
import math
from collections import OrderedDict
from collections import namedtuple
from collections import defaultdict
import itertools
import numpy as np
import pandas as pd
from IPython import embed

from packages.tools_user import *
from frozendict import frozendict


def find_previous_app(app_v, previous_apps):
    """
    caching the supp_index of preivusly supports found by AL for this app

    though this is an importat AL function, I move it to commons to be avaialbe for both the AL module and the MPI module, otherwise i would need the mpi module to see the AL module

    Random must be excluded otherwise it enters in the supervised learning loop
    """
    supp_idx = None
    if previous_apps and '_fg_' not in app_v['supp_type'] and '_rdm_' not in app_v['supp_type']:
        app_previous_AL = [[val.get_prp('supp', 'valid_idx', no_go_2a2=True), val.get_prp(
            'supp', 'size', no_go_2a2=True)] for val in previous_apps.values() if val.do_I_belong(app_v)]
        if app_previous_AL:
            last_app = max(app_previous_AL, key=lambda x: x[1])
            supp_idx = last_app[0]
    return supp_idx


def do_obs_grid(num_dim, over_all, over_some=None, tau_gen=None):
    """
    Generates cartesian observation grid
    """

    # assuming d-cub [0,1]
    if not tau_gen:
        tau_gen = [[0, 1] for _ in range(num_dim)]

    # trimming discretization for non used dimensions
    over_all = [t for t, _ in zip(over_all, range(num_dim))]
    over_some = [t for t, _ in zip(over_some, range(num_dim))]

    # basig obs grid
    ret_grid = [list(np.linspace(ti[0], ti[-1], N)) for ti, N in zip(tau_gen, over_all)]

    # Including data sites tau
    if tau_gen:
        for i, ti in enumerate(tau_gen):
            ret_grid[i] = ret_grid[i] + ti

    # including additional points
    if over_some:
        for i, N in enumerate(over_some):
            min_obs = min(ret_grid[i])
            min_obs = 0
            ret_grid[i] = ret_grid[i] + list(np.linspace(min_obs, 0.1, N))

    # Remove reapeted evaluation points and order the array
    ret_grid = [sorted(set(ti)) for ti in ret_grid]

    return [point for point in itertools.product(*ret_grid)]


def missions_generator(*args, **kwargs):
    """
    Auxiliary function to define the cases to run
    """
    saving_keys = kwargs.pop('saving_keys')
    if kwargs:
        raise ValueError('bad kwargs')

    ret_mission = []
    for mission in [list(mission) for mission in itertools.product(*args)]:
        if not (mission[0] == 'mxs' and mission[1] == 'i-wise'):
            ret_mission.append({key: val for key, val in zip(saving_keys, mission)})

    for mission in ret_mission:
        if mission['sub_kind'] == 'i-wise' or mission['sub_kind'] == 'xs-wise' or mission['kind'] == 'mxs':
            mission.pop('consider_MACRT')

    # eliminate xs-wise for mxs since it doesnt make sense
    ret_mission = [val for val in ret_mission if not (val[
        'kind'] == 'mxs' and val['sub_kind'] == 'xs-wise')]

    user_allowed = []
    for m in ret_mission:
        allow = True
        if m['sub_kind'] == 'i-wise':
            allow = False
        if m['kind'] == 'xs' or m['kind'] == 'mxs':
            # if m['sub_kind']!='all' and m['unit']=='e' :
            #     allow=False

            if m['unit'] == 'e_pcm':
                allow = False
        if m['kind'] == 'k_inf':
            if m['sub_kind'] != 'all':
                allow = False
        if m['kind'] == 'si':
            if m['sub_kind'] != 'all':
                allow = False
            if m['unit'] not in ['e', 'e_per']:
                allow = False
        if allow:
            # print m
            user_allowed.append(m)

    ret_mission = [frozendict(m) for m in user_allowed]  # if package availabe


    return ret_mission



def tensor2line(d, grid, line):
    """
    Projector of tensor. USage:
    prejected_line=[None,0.0,0.2]
    x_value,x_index=tensor2line(dimension,grid,prejected_line)

    change to return x and y value, i.e. add this
        y_g=[]
    for index,point in enumerate(test_dic['test'][common_test][TYPE_DATA][i][r][g]['g(x)'][0]):
            if index in x_index:
                y_g.append(point)
    """

    # the process starts by defining to which dimension reduce the data tensor.
    print d, grid, line
    if d != len(line):
        raise ValueError('of correct type')

    x_value = []
    x_index = []

    d_projection = line
    # I save the point and the index, for latter add more stuff using the index

    for index, point in enumerate(grid):

        # find if the point can be added
        incorporate_flag = 'YES'
        for x_i, d_i in zip(point, d_projection):
            if d_i is None:
                continue
            if d_i != x_i:
                incorporate_flag = 'NO'
                break

        if incorporate_flag == 'YES':

            # notice that a vector is spitted, could be more than one match
            pos = np.where(np.array(d_projection) == None)[0]

            if len(pos) > 1:
                raise ValueError('This should be projection to 1 dimension, not several')
            pos = pos[0]  # from [x] to x

            x_value.append(point[pos])
            x_index.append(index)


    if not x_index:
        raise ValueError("The demanded proyection resulted in an empty x_vec list")
    return x_value, x_index


def do_AX_xs_data(xs_dic, d, aux_tag):
    """
    Read AX data
    """
    AX_dic = A2_read(xs_dic, aux_tag)
    tag = AX_dic['A2'].keys()[0]
    save_obj_pickle(
        AX_dic,
        xs_dic['path']['pickle'],
        'xs_data_' +
        nice(
            AX_dic['A2'].keys()[0],
            'tag') +
        '_d=' +
        str(d) +
        '_' +
        'IRG=' +
        AX_dic['A2'][tag]['info']['IRG'])
    print "Finish extracting A2 data", tag
    # sys.exit()


def order_by(tags, order_list):
    """
    Orders list of string:tags according to the order defined by the
    list of string: order_list
    """
    ret_tag = []
    for order in order_list:
        for tag in tags:
            if all(x in tag for x in order):
                ret_tag.append(tag)


def print_dict(dictionary, ident=''):
    """ Recursively prints nested dictionaries."""

    for key, value in dictionary.iteritems():
        if isinstance(value, dict):
            print ident + '%s' % key
            print_dict(value, ident + '  ')
        else:
            try:
                print_aux = value[0]
                if hasattr(value[0], '__len__'):
                    aux = len(value[0])
                    if aux > 5:
                        print_aux = str(value[0][0:5]).replace(']', '') + '...'
                print ident + '%s = [%s ,...]' % (key, print_aux)
            except IndexError:
                print ident + '%s = %s' % (key, value)


def fgentz(x, f, c, w):
    """
    Dummy function data generator

    d=len(x)=len(c)=len(w)

    input: tuple
    outpu=float
    """

    if len(x) != len(c) != len(w):
        raise ValueError("Fgentz not len(x)!=len(c)!=len(w)")
    if f == 'OSCI':
        # cos in all dimensions
        return np.cos(2 * np.pi * sum([w_i + x_i * c_i for w_i, x_i, c_i in zip(w, x, c)])) + 1
    if f == 'OSCI0':
        # cos in only one dimension for seeing better the grid
        return np.cos(2 * np.pi * (w[0] + x[0] * c[0])) + 1
    if f == 'GAUS':
        return np.exp(-sum([c_i**2 * (x_i - w_i)**2 for w_i, x_i, c_i in zip(w, x, c)]))
    if f == 'CONT':
        return np.exp(-sum([c_i * abs(x_i - w_i) for w_i, x_i, c_i in zip(w, x, c)]))
    if f == 'PPEK':
        return np.product([(c_i**-2 + (x_i - w_i)**2)**-1 for w_i, x_i, c_i in zip(w, x, c)])
    if f == "linear":
        return np.sum([ci * xi for ci, xi in zip(c, x)])
    if f == "linear0":
        return c[0] * x[0]
    raise ValueError("Unknown genz required function")
