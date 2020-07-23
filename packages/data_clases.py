# !/usr/bin/env python
"""
This is a general file for classes to do XS representation

Object-oriented approach
"""
from tools_SG_splines import *
from packages.pppack import S
from packages.tools_SG_common import *
from abc import ABCMeta, abstractmethod
from functools import reduce

import KM.Kfunctions._kernels as KM
import KM.Kfunctions.fortran.fortran_k as FKM
import scipy
from sklearn.model_selection import GridSearchCV
from sklearn.kernel_ridge import KernelRidge
import math as mth
from IPython import embed


# Nested clases can't be pickled, so all name tuples are define as a general class here


class kTupPickable(namedtuple('xs_tup', 'x')):
    """
    Charachterizes OS
    """
    pass

class MxsTupPickable(namedtuple('xs_tup', 'x r g')):
    """
    Charachterizes OS
    """
    pass


class XsTupPickable(namedtuple('xs_tup', 'x i r g')):
    pass


class StateParPickable(namedtuple('grid_point', 'dim_name py_idx a2_idx')):
    pass


class API(object):
    """documentator use as __doc__ = API() in classes"""

    def _print_values(self, obj):
        def _print_value(key):
            if key.startswith('_'):
                return ''
            value = getattr(obj, key)
            if not hasattr(value, 'im_func'):
                doc = type(value).__name__
            else:
                if value.__doc__ is None:
                    doc = 'no docstring'
                else:
                    doc = value.__doc__
            return '%s : %s' % (key, doc)
        res = [_print_value(el) for el in dir(obj)]
        return '\n'.join([el for el in res if el != ''])

    def __get__(self, instance, klass):
        if instance is not None:
            return self._print_values(instance)
        else:
            return self._print_values(klass)


class Approximation(object):
    __metaclass__ = ABCMeta


class Metadata(object):
    """
    I'm putting here a  general class with all path/execution/flag related stuff.
    Class atributes are set with a classmethod.
    """
    @classmethod
    def set_run_date(cls, date):
        """Save the date of the run"""
        cls.run_date = date

    def __init__(self, **kwargs):
        """
        I do need instance atributes for future JdD.
        """
        self.db_I_raw_n = kwargs.pop('db_I_raw_n')
        self.root_input_f = kwargs.pop('root_I_f')
        self.db_f = kwargs.pop('db_f')
        self.root_output_f = kwargs.pop('root_O_f')
        self.num_dim = kwargs.pop('num_dim')
        self.db_O_std_f = kwargs.pop('db_O_std_f')
        self.db_I_raw_f = kwargs.pop('db_I_raw_f')

        self.raw_I_path = Metadata.compute_path(self.root_input_f, self.db_f, self.db_I_raw_f)
        self.std_O_path = Metadata.compute_path(self.root_input_f, self.db_f, self.db_O_std_f)
        for key, val in kwargs.items():
            setattr(self, key, val)

    @staticmethod
    def compute_path(*args):
        """
        Private computation of the path
        """
        return '/'.join(args).replace('///', '/').replace('//', '/')  # use RE

    def __str__(self):
        """
        Return metada info
        """
        return 'str: a2 data coming from: ' + str(self.__dict__)


class MXS_Suite(object):
    """Mix-in class"""
    def variables_mungling1(self):
        """
        this function requires fully developed approximation classes, i.e. with test assign and so on
        here the A2/app dual class use is 'solved'
        Not at all evident that this is a happy solution
        """
        try:#this is a list
            self.mxs_nametuple_space #access from ref
        except AttributeError:
            self.mxs_nametuple_space=self.get_prp('a2','mxs_nametuple_space') #access from app
        try:#this is a list
            self.xs_nametuple_space #access from ref
        except AttributeError:
            #print 'there'
            self.xs_nametuple_space=self.get_prp('a2','xs_nametuple_space') #access from app

    def variables_mungling2(self):
        """
        this function requires fully developed approximation classes, i.e. with test assign and so on
        here the A2/app dual class use is 'solved'
        Not at all evident that this is a happy solution
        """
        try:#this is a function
            self.conc_grid #access from ref
        except AttributeError:
            self.conc_grid=self.get_fprp('test','conc_grid') #access from app
        try:
            self.xs #ref
        except AttributeError:
            def unpickable_f_aux(xs):
                return self.get_eval_dict('xs',xs,unit='e')
            self.xs=unpickable_f_aux

    def compute_mxs(self):
        """
        Computes macroscopic cross sections
        """
        flat_mxs_dict = OrderedDict()
        self.variables_mungling1() # all dirty work to make the function work with a2 class and app class
        self.variables_mungling2() # all dirty work to make the function work with a2 class and app class

        for mxs_name, mxs_tup in self.mxs_nametuple_space.iteritems():

            aux_dict = OrderedDict()
            aux_dict['barn'] = OrderedDict()
            aux_dict['1/cm2'] = OrderedDict()

            for xs in self.find_iso_of_mxs(mxs_tup):
                # this already has the valid index values
                aux_dict['barn'][xs] = np.array(list(self.xs(xs)))  # is a generato
                xs_aux = xs.split('_')[2]
                aux_dict['1/cm2'][xs] = np.array(list(self.conc_grid(xs_aux))) # here the test conc have been hardcoded

            if len(aux_dict['barn'][xs]) != len(aux_dict['1/cm2'][xs]):
                print len(aux_dict['barn'][xs]), len(aux_dict['1/cm2'][xs])
                raise ValueError('mismatch between xs and conc')

            flat_mxs_dict[mxs_name] = MXS_Suite.mxs_calcul(
                aux_dict['barn'], aux_dict['1/cm2'])

        try:
            for xs, val in flat_mxs_dict.iteritems():
                self.set_result('mxs', xs, {'e': val}, type='eval')
        except AttributeError:
            self._flat_mxs_grid_dict = flat_mxs_dict #ApolloContainer does not need to 'incorporate' the data in a specific format
        self.destroy_aux_f()


    def destroy_aux_f(self):
        #cannot be pickled
        if self.xs.__name__=='unpickable_f_aux':
            del(self.xs)
        if self.conc_grid.__name__=='unpickable_f_aux':
            del( self.conc_grid)

    def xs_importance_wrt_mxs(self,key=None):
        """
        Computes the importance of xs
        """
        if not key:
            raise ValueError('what kind of importance do you want?')
        try:
            return self.imp_dict[key]
        except AttributeError:
            def importance(xs_dict, conc_dict, mxs):
                return {xs: np.divide(np.multiply(xs_dict[xs], conc_dict[xs]), mxs) for xs in xs_dict}

            evol_cond = {evol['tuple_idx']: evol['evol_idx'] for evol in self.evol.values()}
            evol_evaluation = lambda x, y: y[x]
            index_grid2evol = [idx for idx, point in enumerate(self.grid_nmlz_idx)
                               if all([evol_evaluation(idx_eval, point) == val_eval for idx_eval, val_eval in evol_cond.iteritems()])]
            if len(index_grid2evol) != len(self.burnup_nmlz):
                raise ValueError('evol vect different than cartesian burnup length')
            proyections={}
            proyections['at_evol'] = lambda x: [x[idx] for idx in index_grid2evol] #proyec to evol cond burnup wise
            proyections['at_grid'] = lambda x: [xi for xi in x] # do nothing, i.e. be at grid position
            self.imp_dict={}
            for proy_n,proy_f in proyections.iteritems():
                #Proyected importance
                imp_dict = defaultdict(OrderedDict)
                for mxs_name, mxs_tup in self.mxs_nametuple_space.iteritems():
                    aux_dict = defaultdict(OrderedDict)
                    for xs in self.find_iso_of_mxs(mxs_tup):
                        aux_dict['barn'][xs] = np.array(proy_f(self.xs(xs)))
                        aux_dict['1/cm2'][xs] = np.array(proy_f(self.conc_grid(xs.split('_')[2])))
                    imp_dict[mxs_name] = importance(
                        aux_dict['barn'], aux_dict['1/cm2'], proy_f(self._flat_mxs_grid_dict[mxs_name]))
                self.imp_dict[proy_n]=imp_dict
            return self.imp_dict[key]

    def default_mxs(self, key):
        """
        defalut get
        """
        return self._flat_mxs_grid_dict[key]

    def find_iso_of_mxs(self, mxs_tup):
        """
        this is nametuple dependent
        this uses mxs in tuple form: xs_tup(x='1', r='abso', g='1')
        """
        mxs_i = []
        # find isotopes than conform the mxs. MACRT excluded
        for xs in self.xs_nametuple_space.values():
            # print xs
            if mxs_tup == tuple([val for val, name in zip(xs, xs._fields)
                             if 'i' not in name]) and xs.i != 'MACRT':
                mxs_i.append(xs)
        # need xs as str and the corresponding isotope
        return [MXS_Suite.build_name(xs) for xs in mxs_i]

    @staticmethod
    def build_name(xs):
        return ''.join(['_' + el for el in xs])

    @staticmethod
    def mxs_calcul(xs_dict, conc_dict):
        """
        The condition 'if dictxs_xs_name != dictconc_xs_name:'
        Is VERY restrictive as it forces the two dict to be order dicts whith the same keys and so one. More general things could be written with proper tests.
        """

        macro = 0
        # numpy operations
        for dictxs_xs_name, dictxs_xs_val, dictconc_xs_name, dictconc_xs_val in zip(
                xs_dict.keys(), xs_dict.values(), conc_dict.keys(), conc_dict.values()):
            if dictxs_xs_name != dictconc_xs_name:
                print dictxs_xs_name
                print dictconc_xs_name
                raise ValueError('Conc does not correspond to xs')
            # broadcast and ele-wise mult
            macro = macro + np.multiply(dictxs_xs_val, dictconc_xs_val)
        return macro


class K_Suite(object):
    """
    Suite of infinite multiplication factor
    """
    def compute_k(self,kind,and_ro=False):
        """
        a proper use of self.mxs_nametuple_space.iteritems() would permit to do this region by region
        Managment by region done by hand. If actual region homogenization is present this will need to be re-written
        this is writting for the app
        """
        aux_mxs_dict={}
        for mxs in ['_1_abso_1', '_1_tran012_1', '_1_nufi_1', '_1_tran021_2', '_1_nufi_2', '_1_abso_2']:
            if mxs not in self.get_prp('a2','mxs_namespace'):
                raise ValueError('maybe just return None somewhere. In any case k calcul not possible')
            aux_mxs_dict[mxs]=self.get_eval_dict('mxs',mxs,unit='e')
        results_dict=self.k_calcul(kind,aux_mxs_dict,and_ro)
        def aux_give_name(st):
            if 'k_inf' in st:
                return 'k_inf'
            if 'si' in st:
                return 'si'
            raise ValueError('strange kind')
        [self.set_result(aux_give_name(kind_i), kind_i, {'e': val}, type='eval') for kind_i,val in results_dict.iteritems()]
    def compute_k_inf_rb(self,kind,and_ro=False):
        """
        this is writting for a2 to avoid using its k_inf
        #!# this and the upper one could and should be fusioned
        """
        aux_mxs_dict={}
        for mxs in ['_1_abso_1', '_1_tran012_1', '_1_nufi_1', '_1_tran021_2', '_1_nufi_2', '_1_abso_2']:
            if mxs not in self.mxs_namespace:
                raise ValueError('maybe just return None somewhere. In any case k calcul not possible')
            aux_mxs_dict[mxs]=self.mxs(mxs,unit='e')
        self.k_inf_rb=self.k_calcul(kind,aux_mxs_dict,and_ro)
    @staticmethod
    def k_calcul(kind,mxs_dict,and_ro=False):
        if kind=='k_inf_classic':
            # the autovector is the SI, a menos de una constante.
            SI_vec = [(sg_a2 + sg_21) / sg_12 for sg_a2, sg_21,
                      sg_12 in zip(mxs_dict['_1_abso_2'],mxs_dict[ '_1_tran021_2'],mxs_dict['_1_tran012_1'])]
            k_vec = [(nufi1 * SI + nufi2) / ((sg_a1 + sg_12) * SI - sg_21) for nufi1, nufi2, SI, sg_12, sg_21, sg_a1 in
                     zip(mxs_dict['_1_nufi_1'], mxs_dict['_1_nufi_2'],SI_vec, mxs_dict['_1_tran012_1'], mxs_dict['_1_tran021_2'], mxs_dict['_1_abso_1'])]
            ro_vec = []
            return {'_1_k_inf_rb': np.array(k_vec,order='F'),'_1_si_rb': np.array(SI_vec,order='F')}
        raise ValueError('unknow kind of k calculation')

    @staticmethod
    def build_name(x,k):
        return ''.join(['_' + el for el in x])+k+'_rb'

class ApolloContainer(MXS_Suite,K_Suite):
    """
    Class of A2 data container and manipulator. Each instance corresponds to a use case. Inner data representation in dic only
    Class is read-only.

    Responsabilities:
    * Contain and exhibit standarized a2 data. Mainly xs  data and the state space in several ways (cartesian, grid, indexes, etc)
    * Saving in database
    * Loading of database.
    """
    identy = 0

    def __init__(self, a2_data, num_dim):
        """Initialization of container"""
        self._a2_data = a2_data
        self._num_dim = num_dim
        self._train_cartesian_idx = None

        # the evol vector should be given, here is estimated
        evol = OrderedDict()
        for dim_tup, dim_val in a2_data['PS']['real']['data'].iteritems():
            if dim_tup.dim_name != 'BURNUP':
                evol[dim_tup.dim_name] = {}
                evol[dim_tup.dim_name]['evol_idx'] = range(len(dim_val))[int(len(dim_val) / 2)]
                evol[dim_tup.dim_name]['tuple_idx'] = dim_tup. py_idx
        self.evol = evol

        # checking format compliance in a few random names.
        for name in [a2_data['xs_names'][i]
                     for i in [random.randint(0, len(a2_data['xs_names']) - 1) for _ in xrange(5)]]:
            if len(name.split('_')) != 5:
                raise RuntimeError('Better format required at standarization')
        self.XIRG = a2_data['XIRG']  # !# output space, any notion of order should be here

        # mxs space generated from available mxs
        mxs_output_space = []
        for xs in self.xs_nametuple_space.values():
            mxs_output_space.append(MxsTupPickable(xs.x, xs.r, xs.g))

        mxs_output_space = set(mxs_output_space)
        mxs_names = [MXS_Suite.build_name(mxs) for mxs in mxs_output_space]
        self._mxs_namespace = mxs_names
        self._mxs_nametuple_space = {name: tup for name, tup in zip(mxs_names, mxs_output_space)}

        k_output_space=[]
        for tup in self.mxs_nametuple_space.values():
            k_output_space.append(kTupPickable(tup.x))
        k_output_space=set(k_output_space)

        k_i='k_inf'
        k_inf_names=[K_Suite.build_name(k,'_'+k_i) for k in k_output_space]
        self.k_inf_namespace=k_inf_names
        self.k_inf_nametuple_space={k_i:tup for k_i,tup in zip(k_inf_names,k_output_space)}
        k_i='si'
        si_names=[K_Suite.build_name(si,'_'+k_i) for si in k_output_space]
        self.si_namespace=si_names
        self.si_nametuple_space={k_i:tup for k_i,tup in zip(si_names,k_output_space)}


    def spezializate_run_name(self):
        ret_val = ['_eff']
        for key, sub_dict in self._a2_data['XIRG'].iteritems():
            ret_val.append(key + str(sub_dict['N']))
        return''.join(ret_val)


    def compute_train_gen_bybudget(self, serma=True, training_budget=None, distri='uniform', less_sermapoints=None):
        """
        Computation of cartesian training generator. The proces is index based. Only Cartesian construction
         are used in the training, the test is the complement of such space, totalizing the
          apollo2 data
        """

        # include border in Cartesian generator of training set
        _train_cartesian_idx = [[tau_i[0], tau_i[-1]]
                                    for tau_i in self._a2_data['PS']['nmlz']['indx']['tau']]

        unit=next(iter(training_budget)) #fastest way to find single key dict
        if unit=='per':
            training_budget_N=[int(len(v)*bg_per/100) for v,bg_per in zip(self._a2_data['PS']['nmlz']['indx']['tau'],training_budget[unit])]
        if unit=='len(idx)':
            training_budget_N=training_budget[unit]

        #serma=False
        # include serma points if required
        if serma:
            self._serma_cartesian_idx = [closest_serma_disct(dim, vec, less_sermapoints) for dim, vec in self._a2_data[
                'PS']['real']['data'].iteritems()]  # finding serma points
            _train_cartesian_idx = [sorted(set(training_i + serma_i)) for
                                        training_i, serma_i in zip(_train_cartesian_idx,
                                                                   self._serma_cartesian_idx)]  # adding serma points
            for i, budget in enumerate(training_budget_N):
                while len(_train_cartesian_idx[i]) >= budget:
                    print 'serma points > training, are you sure?',i,len(_train_cartesian_idx[i]),budget
                    middle = int(len(_train_cartesian_idx[i]) / 2.0)
                    _train_cartesian_idx[i] = _train_cartesian_idx[i][
                        0:middle] + _train_cartesian_idx[i][middle + 1:]

        maybe_unordered_dimnames=[[val,val.py_idx] for val in self._a2_data['PS']['real']['data'].keys()]
        ordered_dimnames=[var[0] for var in sorted(maybe_unordered_dimnames,key=lambda x:x[1])]
        # exhaust available training_budget_N if required
        for training, cartesian_di_idx,dimname, budget in zip(_train_cartesian_idx,
                                          self._a2_data['PS']['nmlz']['indx']['tau'],ordered_dimnames, training_budget_N):
            size_novo = budget - len(training)
            if size_novo < 0:
                raise RuntimeError('Avialiable budget-serma points=%s ' % size_novo)
            random.seed(0)
            #print training, data, budget,size_novo
            #sys.exit()
            while size_novo > 0:
                # the complement of the intersection between data and training is found
                available_data_idx = [point for point in cartesian_di_idx if point not in training]
                index_novo = random.randint(0, len(available_data_idx) - 1) #random choosing
                training.append(available_data_idx[index_novo])
                size_novo -= 1

        self._train_cartesian_idx = [sorted(set(t)) for t in _train_cartesian_idx]

    @property
    def train_cartesian_idx(self):
        return self._train_cartesian_idx

    def train_cartesian_grid_val(self,xs_n=None):
        xs=None
        if xs_n:
            xs=[self.xs(xs_n)[self.grid_nmlz_idx.index(train_idx)] for train_idx in itertools.product(*self._train_cartesian_idx)]
        return [self.grid_nmlz_val[self.grid_nmlz_idx.index(train_idx)] for train_idx in itertools.product(*self._train_cartesian_idx)],xs

    @property
    def train_grid_idx(self):
        return [self.grid_nmlz_idx.index(train_idx) for train_idx in itertools.product(*self._train_cartesian_idx)]

    @property
    def serma_cartesian_idx(self):
        return self._serma_cartesian_idx

    def compute_supp_gen(self, supports_budget):
        """
        Computation of supports from training points for supports_budget. Borders of thetraining are always included
        If possible, serma points are selected first. Big support include smaller ones.
         Only computation of cartesian generation is required
        """
        #supports_budget=cp.deepcopy(supports_budget)
        serma_flag=supports_budget.pop('serma_supp')
        if len(supports_budget)>1:
            raise ValueError('too many notes')
        unit=next(iter(supports_budget)) #fastest way to find single key dict
        if unit=='per':
            supports_budget_N=[]
            for supp in supports_budget[unit]:
                supports_budget_N.append([int(len(v)*bg_per/100) for v,bg_per in zip(self._train_cartesian_idx,supp)])
        if unit=='len(idx)':
            supports_budget_N=supports_budget[unit]

        # chaking that at least there is place for the borders
        for supp in supports_budget_N:
            for supp_i in supp:
                if supp_i<2:
                    raise ValueError("will not be able to contain both borders %s", str(supp))

        # cheking that training points actually permit those percentage without repetition
        if len(supports_budget_N)!=len(set(map(tuple, supports_budget_N))):
            raise ValueError('one of the supports was eliminated %s',str(supports_budget_N))

        # budget must be in increasing order for every dimension
        for budget in zip(*supports_budget_N):
            if sorted(budget) != list(budget):
                raise ValueError('Non-ascending budgets %s' % supports_budget_N)

        # budgent must be coherent with number of dimension
        for budget in supports_budget_N:
            if len(budget) != self._num_dim:
                raise ValueError('dimensional mismatch between budget and number of dimensions')

        # computing support generators
        _support_generators = []
        for budget in supports_budget_N:
            supp = []
            for budget_dim, dim_idx, training_dim, serma_dim in zip(budget, range(len(budget)),
                                                                    self._train_cartesian_idx,
                                                                    self._serma_cartesian_idx):
                size_novo = cp.copy(budget_dim)
                random.seed(0)
                supp_dim = []

                # Adding point until exhaustion of budget
                while size_novo > 0:
                    # trying to use serma points
                    available_data = list(filter(lambda x: x not in supp_dim, serma_dim))

                    if not available_data:
                        try:  # trying to use previous support [-1] for the dimension dim_idx
                            available_data = list(
                                filter(lambda x: x not in supp_dim, _support_generators[-1][dim_idx]))
                        except IndexError:
                            pass

                        if not available_data:
                            # trying to use remaining training points
                            available_data = list(filter(lambda x: x not in supp_dim, training_dim))

                    if not available_data:
                        raise RuntimeError('all serma, and training point have being use,\
                                   but still requiring % s' % size_novo)
                    index_novo = random.randint(0, len(available_data) - 1)

                    # Enforcing inclusion of interval border of the training set in candidate index_novo
                    if training_dim[0] not in supp_dim:  # first element
                        index_novo = available_data.index(training_dim[0])
                    if training_dim[-1] not in supp_dim:
                        index_novo = available_data.index(training_dim[-1])

                    supp_dim.append(available_data[index_novo])
                    size_novo -= 1

                if len(set(supp_dim)) != len(supp_dim):
                    raise RuntimeError('Redundancy in construction of the support space')

                if training_dim[0] not in supp_dim or  training_dim[-1] not in supp_dim:
                    raise RuntimeError('the budget was so small that could not include the borders %s', str(budget))

                supp.append(sorted(supp_dim))

            _support_generators.append(supp)

        if serma_flag=='Y':
            _support_generators.append(self._serma_cartesian_idx)

        self._support_generators = _support_generators

    @property
    def supp_gen_idx(self):
        """
        The supports are generated at this Apollo level
        """
        return self._support_generators

    def centralized_dict_access(self, key):
        """Its better to have a few centralize chooosings like this or
        each propery defines its access?"""
        if key == 'burnup_rel':
            for name_tup in self._a2_data['PS']['real']['data']:
                if name_tup.dim_name == 'BURNUP':
                    break

            return self._a2_data['PS']['real']['data'][name_tup]

        if key == 'burnup_nmlz':  # !# this is cart! but it is not indicated
            for name_tup in self._a2_data['PS']['nmlz']['tau_gen']:
                if name_tup.dim_name == 'BURNUP':
                    break

            return self._a2_data['PS']['nmlz']['tau_gen'][name_tup]

    @property
    def burnup_real(self):
        key = 'burnup_real'
        return self.centralized_dict_access(key)

    def conc_cartesian(self, key):
        # this could be the standard way to acess any usefull piece of data
        if key in self._a2_data['PS']['real']['conc']:
            return self._a2_data['PS']['real']['conc'][key]
        print self._a2_data['PS']['real']['conc'].keys()
        raise ValueError('key %s not found in conc dict', key)

    def conc_grid(self, key):
        # this could be the standard way to acess any usefull piece of data
        # in principle the key='U235' or the isotope directly
        #!# note that they are not stored by region but only by isotope
        if key not in self._a2_data['PS']['real']['conc']:
            dynamic_rekeying=[part in self._a2_data['PS']['real']['conc'] for part in key.split('_')] # to also (try to) catch _1_AM241_abso_1
            dynamic_rekeying=[idx for idx,val in enumerate(dynamic_rekeying) if val==True]
            if len(dynamic_rekeying)>1:
                raise ValueError('dynamic rekying got messed up')
            key_idx=dynamic_rekeying[0]
            key=key.split('_')[key_idx]
            if key not in self._a2_data['PS']['real']['conc']:
                raise ValueError('key %s not found in conc dict after dynamic rekeyng', key)
        try:
            self.bu_idx
        except AttributeError:
            self.bu_idx = [
                dim.py_idx for dim in self.tuple_index_order if dim.dim_name == 'BURNUP'][0]

        return [self._a2_data['PS']['real']['conc'][key][idx[self.bu_idx]] for idx in self.grid_nmlz_idx]

    @burnup_real.setter  # does this has a real utility?
    def burnup_real(self, value):
        self.burnup_real = self.burnup_real

    @property
    def burnup_nmlz(self):
        key = 'burnup_nmlz'
        return self.centralized_dict_access(key)

    @burnup_nmlz.setter  # does this has a real utility?
    def burnup_nmlz(self, value):
        self.burnup_nmlz = self.burnup_nmlz

    def xs(self, name):
        """dictionary based XS access"""
        if len(name.split('_')) == 5:
            _, x, i, r, g = name.split('_')
            try:
                return self._a2_data['xs'][x][i][r][g]
            except KeyError:
                raise AttributeError('requested isotope not found: %s' % name)
        # this being an AtributeError is vital, as is the expected value of a non found attr
        # by __setstate__ in the unpickling process
        raise AttributeError('the attributed "%s" is not a a xs xirg' % name)

    def mxs(self,name,unit='e'):
        try:
            if name in self._flat_mxs_grid_dict:
                return self._flat_mxs_grid_dict[name]
        except AttributeError:
            raise ValueError('mxs in dict or mxs_dict not present')


    def k_inf(self,name,unit='k_inf'):
        if name=='_1_k_inf_a2':
            """
            This is the a2 data set, it can be somewhat different than what calculated with formulas
            """
            return self._a2_data['k']['k_inf'][unit]
        if name=='_1_k_inf_rb':
            """
            This uses the formular of k_inf utilizing rebuild MXS only. If this coinsided with than formulation using MACRT depends on the proper build of the MXS
            and the considering of all the isotopes and solving the RES nufi problem. Its safer for the analysis
            If nufi RES is not consider this will differenciate with its calculation with MACRT in about 100pcm and in the a2 data in about 500 pcm for unknown reasons
            """
            try:
                return self.k_inf_rb['_1_k_inf_rb']
            except AttributeError:
                raise ValueError('problems with k rebuit')


    def si(self,name):
        """
        Not that this explicitly jumps over the whole x definition, it will need to be re-writting if x actually changes
        """
        if name=='_1_si_a2':
            return self._a2_data['fi']['si']['si']
        if name=='_1_si_rb':
            try:
                return self.k_inf_rb['_1_si_rb']
            except AttributeError:
                raise ValueError('problems with k rebuit')

    @property
    def a2_ps_cardinality(self):
        return self._a2_data['PS']['nmlz']['N_tau']

    @a2_ps_cardinality.setter
    def a2_ps_cardinality(self, value):
        raise RuntimeError('a2 data read-only')

    @property
    def ref_2_a2(self):
        return self._a2_data

    @ref_2_a2.setter
    def ref_2_a2(self, value):
        raise RuntimeError('a2 data read-only')

    @property
    def tuple_index_order(self):
        return self._a2_data['PS']['real']['data'].keys()

    @property
    def input_space_names(self):
        return [val.dim_name.lower() for val in self._a2_data['PS']['real']['data'].keys()]

    @property
    def xs_namespace(self):
        return self._a2_data['xs_names']

    @xs_namespace.setter
    def xs_namespace(self, value):
        raise RuntimeError('a2 data read-only')

    @property
    def xs_nametuple_space(self):
        return self._a2_data['xs_tuple']

    @xs_nametuple_space.setter
    def xs_nametuple_space(self, value):
        raise RuntimeError('a2 data read-only')

    @property
    def mxs_namespace(self):
        return self._mxs_namespace

    @property
    def mxs_nametuple_space(self):
        return self._mxs_nametuple_space

    @property
    def grid_nmlz_val(self):
        return self._a2_data['PS']['nmlz']['grid']

    @grid_nmlz_val.setter
    def grid_nmlz_val(self, value):
        raise RuntimeError('a2 data read-only')

    @property
    def grid_nmlz_idx(self): #!# grid, tau, gen etc is not a happy notation. it should be A2 values there are many levels, A2,TT,sup, etc and this brings confusion
        """
        #!# if it has idx then is the same if normalized or not, idx is another unit: idx, nmlz, real
        """
        return self._a2_data['PS']['nmlz']['indx']['grid']

    @grid_nmlz_idx.setter
    def grid_nmlz_idx(self, value):
        raise RuntimeError('a2 data read-only')

    @property
    #tau_gen is so unclear usea apollo2cartesian grid or something like that
    def tau_gen_val(self):
        return self._a2_data['PS']['nmlz']['tau_gen']

    @tau_gen_val.setter
    def tau_gen_val(self, value):
        raise RuntimeError('a2 data read-only')

    def __str__(self):
        return '(str) a2 container: ' + str(self._a2_data['XIRG'])

    def __repr__(self):
        return '(repr) a2 container:' + str(self._a2_data['XIRG'])


class Meta_Methods(object):
    """
    Methods present present in Mediator and Approximation
    """
    pass

class Mediator(Meta_Methods):
    """Support class. Hopefully I manege to only hodl a reference to the memory location of
     the actual points. It only needs to see some part of the instance a2_data and store the correct references
     Valid index method only usable for vectorized data.

    Responsabilitie:
    * Be intermediary between data and approximation, without auxiliary memory allocations
    * gives supp points

    note that this class is cartesian-free can be used for a scattered suppoort
     """

    def __init__(self, a2_inst):
        self._ref_to_data = a2_inst  # this should be only a reference
        self._valid_idx = None


    def Cartesian_2_iterator_core(self):
        """ A list of valid indexes for this support is build"""
        # finding valid vectorization index for support's cartesian construction
        valid_idx = []
        #valid_grid_idx = []
        if not self.cartesian_supp_index:
            raise RuntimeError('Cannot compute an interator core without a cartesian training set')

        for index, grid_point in enumerate(self._ref_to_data.grid_nmlz_idx):
            belong_to_training = [x_dim in cartesian_supp_index_dim for x_dim,
                                  cartesian_supp_index_dim in zip(grid_point, self.cartesian_supp_index)]
            if self.boolean_TTcondition(belong_to_training):
                valid_idx.append(index)
        self._valid_idx = valid_idx
        self.compute_bu_valid_idx()


    def compute_bu_valid_idx(self): 
        bu_idx = [dim.py_idx for dim in self.get_prp_valid('tuple_index_order',v_bypass=True) if dim.dim_name == 'BURNUP'][0]
        self._valid_grid_idx=[self._ref_to_data.grid_nmlz_idx[valid] for valid in self._valid_idx]
        self.bu_valid_idx = sorted(
            set([el[bu_idx] for el in self._valid_grid_idx]))  # forploting conc


    def valid_iterator(self, data): 
        """ Generates a valid iterator.
        In theory the data is only a reference without duplication and
         yield produces stream as needed. Only vectorized data
        """
        for valid_index in self._valid_idx:
            yield data[valid_index]

    @property
    def im_serma(self):
        if all([d=='S' for d in self._im_serma]):
            return True
        else:
            return False

    def get_prp_valid(self,kind,v_bypass=False,no_go_2a2=False,sub_kind=False):
        """
        This returns the property passed by the valid
        """
        if no_go_2a2: # to access properties of the class itself
            return getattr(self,kind)
        values=getattr(self._ref_to_data,kind)
        if v_bypass!=False:
            return values
        else:
            if 'multiple_xs_valid_iterators' in self.__dict__.keys() and 'grid' in kind:
                return self.dict_valid_iterator(values,sub_kind)
            else:
                #normal case
                return self.valid_iterator(values)

    def get_fprp_valid(self,mthd,kind=False,sub_kind=False):
        """
        This returns the property of the ref_to_data, comming from a FUNCTION and then applying valid.
        If no kind/subkind if given then it is superinteded that the FUNCTION itself is wanted and so is given
        """
        if kind!=False:
            method=getattr(self._ref_to_data,mthd)
            if sub_kind!=False:
                values=method(kind,sub_kind)
            else:
                values=method(kind)

            if 'multiple_xs_valid_iterators' in self.__dict__.keys() and mthd=='xs':
                return self.dict_valid_iterator(values,kind)
            else:
                #normal case
                return self.valid_iterator(values)
        else:
            # get_fprp_valid returns f(arg) if no arg is given, then is intended that 'f' is wanted
            def unpickable_f_aux(kind,sub_kind=False):
                """
                When returning the original function is necessary to provide the valid iterator as well
                """
                method=getattr(self._ref_to_data,mthd)
                return self.valid_iterator(method(kind))
            return unpickable_f_aux


    def conc_eval_cart(self, name):
        return [self._ref_to_data.conc_grid(name)[val] for val in self.bu_valid_idx]

    @property
    def valid_grid_idx(self):
        """
        [(1, 0), (3, 0), (4, 0), (5, 0), (7, 0), (9, 0), (10, 0), (18, 0)...
        """
        return self._valid_grid_idx

    @property
    def valid_idx(self):
        return self._valid_idx

    @property
    def burnup_nmlz_cart(self):
        """
        gives back the burnup vector seen by the support. IT provides an alternative valid iterator obtanied from the grid
        """
        return [self._ref_to_data.burnup_nmlz[val] for val in self.bu_valid_idx]

    def __str__(self):
        return str(self._ref_to_data)

    @property
    def cartesian_val(self):
        return self._cartesian_val

    @property
    def cartesian_supp_index(self):
        return self._cartesian_supp_index

    @cartesian_supp_index.setter  
    def cartesian_supp_index(self, value):
        self._cartesian_supp_index = value


class Test(Mediator): 
    """
    Responsabilitie:
    * Support instanziation and defining the boolean condition for belonguing to the set
    """

    def __init__(self,a2_inst, cartesian_training_index=None, valid_idx=None):
        """
        cartesian_supp_index is the training set. Name should be change acordanly and any conflict in mediator solved #!#
        """
        super(Test, self).__init__(a2_inst)
        if cartesian_training_index:
            self.cartesian_supp_index=cartesian_training_index

        if valid_idx:
            self._size=len(valid_idx)
            self._valid_idx=valid_idx
            self.compute_bu_valid_idx()

        if not cartesian_training_index and not valid_idx:
            raise RuntimeError("cant instanciate Test")

    @property
    def size(self):
        """
        The test has a clear definition: amout of test points
        """
        return self._size


    @staticmethod
    def boolean_TTcondition(belong_to_training):
        """if any  element of the grid point does not belong to the cartesian support, is a valid test point"""
        return any([not b for b in belong_to_training])

    def __str__(self):
        return 'Im test of size= '+str(len(self._valid_idx))

    def sparsefy_iterator_core(self,density_regulizer,test_lim=False,bound_to_training=False,triche_grid='NO'):
        """
        User to eliminate elements of valid iterator core.
        This is a simplistic method, a proper histrogram  normalization should be done
        Though the valid iterator core is modiffied, the test is no longuer cartesian
        """

        bu_lim=0.015 # here is the biggest concentration
        i=0
        logic_flag=True
        # eliminates any training point that is not in a line where training points are avaialbe. This avoids offset errors
        if bound_to_training!=False:
            Lold=len(list(self._valid_idx))
            #for the indeces that are not burnup
            non_burnup_idx_py=[tup.py_idx for tup in self.get_prp_valid('tuple_index_order',v_bypass=True) if tup.dim_name!='BURNUP']
            _valid_idx_new=[]
            for grid_idx,valid_idx in zip(self.get_prp_valid('grid_nmlz_idx'),self._valid_idx):
                #check that every index of a grid point(except burnup) belongs to the cartesian grid
                if all([grid_idx[i] in list(self.cartesian_supp_index)[i] for i in non_burnup_idx_py]):
                    _valid_idx_new.append(valid_idx)
            self._valid_idx=_valid_idx_new
            print 'by bound to treining reduction in test points from ', Lold, 'to ',len(_valid_idx_new)


        if density_regulizer>0:
            Lold=len(list(self._valid_idx))
            # this reduced the amount of points if Bu<0.015 by doing one YES one NO
            while i<density_regulizer:
                _valid_idx_new=[]
                for point,valid_idx in zip(self.get_prp_valid('grid_nmlz_val'),self._valid_idx):
                    if point[0]>bu_lim:
                        _valid_idx_new.append(valid_idx)
                        continue
                    if point[0]<bu_lim and logic_flag==True:
                        _valid_idx_new.append(valid_idx)
                        logic_flag=False
                    else:
                        logic_flag=True
                i=i+1
                self._valid_idx=_valid_idx_new
            print 'by density regulizer in test points from ', Lold, 'to ',len(_valid_idx_new)

        if triche_grid=='YES':
            Lold=len(list(self._valid_idx))
            #for the indeces that are not burnup
            bu_idx_py=[tup.py_idx for tup in self.get_prp_valid('tuple_index_order',v_bypass=True) if tup.dim_name=='BURNUP'][0]
            _valid_idx_new=[]
            for grid_val,valid_idx in zip(self.get_prp_valid('grid_nmlz_val'),self._valid_idx):
                #check that every index of a grid point(except burnup) belongs to the cartesian grid
                if grid_val[bu_idx_py]>0.2:
                    _valid_idx_new.append(valid_idx)
            self._valid_idx=_valid_idx_new
            print 'by triche grid reduction in test points from ', Lold, 'to ',len(_valid_idx_new)

        # I consider that more than 10000 test points is not necessary so
        if test_lim!=False:
            Lold=len(list(self._valid_idx))
            _valid_idx_new=cp.deepcopy(self._valid_idx)
            logic_flag=True
            while len(_valid_idx_new)>test_lim:
                for i,_ in enumerate(_valid_idx_new):
                    if logic_flag==True:
                        del _valid_idx_new[i]
                        logic_flag=False
                    else:
                        logic_flag=True
                    i=i+1
            self._valid_idx=_valid_idx_new
            print 'by test lim in test points from ', Lold, 'to ',len(_valid_idx_new)

        if len(self._valid_idx)!=len(set(self._valid_idx)):
            raise ValueError('A duplication of test points ocurred in test grid manipulation')
        if len(self._valid_idx)!=len(list(self.get_prp_valid('grid_nmlz_val'))):
            raise ValueError('broken test class: mismatch between valid index and grid')

        self.compute_bu_valid_idx()


class Support(object):
    @property
    def size(self):
        """
        in the support the isze is the amount of coefficients
        """
        if 'multiple_xs_valid_iterators' in self.__dict__.keys():
            size=sum(len(coef) for coef in  self._valid_xsw_idx.values())
        else:
            size=len(self._valid_idx)*len(self.get_prp_valid('xs_namespace',v_bypass='true'))
        return size

class Scattered_support(Mediator,Support):
    """
    does not have a cartesian notion. It recived directly the valid index
    """
    def __init__(self, a2_inst,name,valid_idx=None):
        super(Scattered_support, self).__init__(a2_inst)
        self.name=name
        #!# since the use of __getattr__ sends me directly to the a2container, I need to explicitly say that this is None here. Only for this reason that defualt behaivour should be removed
        #!# then it would be possible use try excelt Atribute error in Mediator properties
        self._cartesian_val=None #boilerplate code
        self._cartesian_supp_index=None #boilerplate code
        self._im_serma=[False]

        if valid_idx:
            if type(valid_idx)==list:
                self._valid_idx=valid_idx
                self.compute_bu_valid_idx()

            if type(valid_idx)==dict:
                self._valid_xsw_idx=valid_idx
                self.multiple_xs_valid_iterators=True

    def dict_valid_iterator(self, data,xs): #,other_valid_idx=None
        """ Generates a valid iterator.
        In theory the data is only a reference without duplication and
         yield produces stream as needed. Only vectorized data
        """
        for valid_index in self._valid_xsw_idx[xs]:
            yield data[valid_index]


    def __str__(self):
        return 'Scattered support of N='+str(self.size)

class Cartesian_support(Mediator,Support):
    """
    Responsabilitie:
    * Support instanziation and defining the boolean condition for belonguing to the set
    * Exposing cartesian values and indexes

    unit test: that length of requested valid points coicides with length supoort
    This class supports cartesian logic
     """
    __doc__ = API()

    def __init__(self, cartesian_supp_index, a2_inst,name,cartesian_flag=True):
        super(Cartesian_support, self).__init__(a2_inst)
        self._cartesian_supp_index = cartesian_supp_index
        self.name = name

        # define which dimensions are denser than serma disct
        self._im_serma = []
        for supp_dim_gen, serma_dim_gen in zip(
                cartesian_supp_index, self._ref_to_data.serma_cartesian_idx):
            serma_points_presence=all(x in supp_dim_gen for x in serma_dim_gen)
            if serma_points_presence:
                if len(supp_dim_gen)==len(serma_dim_gen):
                    self._im_serma.append('S')
                else:
                    self._im_serma.append('B')
            else:
                self._im_serma.append('F')

        self._cartesian_val = []
        for cartesian_supp_index_dim, tau_gen_val_dim in zip(
                cartesian_supp_index, a2_inst.tau_gen_val.values()):
            self._cartesian_val.append([tau_gen_val_dim[i] for i in cartesian_supp_index_dim])

    @property
    def budget(self):
        return len(list(itertools.product(*self.cartesian_supp_index)))

    @staticmethod
    def boolean_TTcondition(belong_to_training):
        """if the grid point belongs to the cartesian support for all dimensions, is valid index"""
        return all(belong_to_training)

    def __str__(self):
        return 'str Cartesian_support ' + self.name

    def __repr__(self):
        return 'repr Cartesian_support ' + self.name

class Statistics_Metaclass(object):
    pass

class Approximation_Container(Statistics_Metaclass, MXS_Suite,K_Suite, Meta_Methods):
    """
    This class contains the approximation and evaluation and statistical data
    The evaluation and statistics are done outside the class sicne are stateless.

    content:
    *evaluated xs
    *approximation's errors
    *statistical result of approximation's errors
    """
    __metaclass__ = ABCMeta
    __doc__ = API()

    def __init__(self, supp,app_n,supp_type):
        self._supp_inst = supp  # here this should be only a reference
        self._test_inst = None  # not very sure. But i need the conc at test, and to preserve the interface and encapsulation I need the instance reference
        self._xirg = supp.get_prp_valid('XIRG',v_bypass=True)

        self.dict_common_app={}
        self.dict_common_app['supp_type']=supp_type
        self.dict_common_app['size']=supp.size

        self._evaluation_dict = OrderedDict()
        self._obs_dict = OrderedDict()
        self._error_dict = OrderedDict()
        self._statistics_dict = OrderedDict()

        self._present_kinds = ['xs', 'mxs','k_inf','si'] #!# automatize this
        _names = OrderedDict()
        for kind in self._present_kinds:
            tup = self._supp_inst.get_prp_valid(kind + '_nametuple_space',v_bypass=True).values()[0]
            _names[kind] = [el.replace('region', 'x').replace('iso', 'i').replace(
            'reac', 'r').replace('group', 'g') for el in list(tup._fields)]
        self.OS_kind_names = _names
        self._xs_os_order = [os.lower() for os in self._xirg]
        self.OS_kind_names

    def __str__(self):
        return '(str) app container: ' + str(self._present_kinds)


    @property
    def im_serma(self):
        return self._supp_inst.im_serma

    def get_prp(self,who,kind,no_go_2a2=False,sub_kind=False):
        """
        accessing the supp or test structures
        """
        if who=='supp':
            return self._supp_inst.get_prp_valid(kind,no_go_2a2=no_go_2a2,sub_kind=sub_kind)
        if who=='test':
            return self._test_inst.get_prp_valid(kind,no_go_2a2=no_go_2a2,sub_kind=sub_kind)
        if who=='a2':
            return self._supp_inst.get_prp_valid(kind,v_bypass=True)
        raise ValueError('dont know where u want the data from')

    def get_fprp(self,who,mthd,kind=False,sub_kind=False):
        """
        here only supp or test as this values pass through valid. For raw information as the number of I or XIRG use property/who='a2'
        """
        if who=='supp':
            h_data=self._supp_inst
        if who=='test':
            h_data=self._test_inst

        return h_data.get_fprp_valid(mthd,kind=kind,sub_kind=sub_kind)
        print who,mthd,kind,sub_kind
        raise ValueError('dont know where u want the data from')

    def conc_supp_cart(self, key):
        return self._supp_inst.conc_eval_cart(key)

    def get_eval_dict(self,kind,sub_kind,unit='e'):
        """
        accessing the eval dict
        """
        return self._evaluation_dict[kind][unit][sub_kind]

    @property
    def types(self):
        if all(t == self._types[0] for t in self._types):
            return 'all_' + self._types[0]
        return self._types

    @types.setter
    def types(self, value):
        self._types = value

    @property
    def orders(self): #!# very ambiguos with order of tuple in inputspace
        if all(t == self._orders[0] for t in self._orders):
            return 'all_' + str(self._orders[0])
        return str(self._orders)

    @orders.setter
    def orders(self, value):
        self._orders = value


    def str_2_idy(self, xs_name):
        """
        xs in general are flattened here, and for fast querring the corresponding
        name tuples are available for each xs
        """
        _, x, i, r, g = xs_name.split('_')
        return XsTupPickable(x, i, r, g)

    def get_xs_id(self, key):
        return self._approx_dic[key]['xs_id']

    def get_name(self,key):
        if key=='math':
            return self.dict_special_app['math_name']
        if key=='app':
            return self.dict_special_app['app_name']
        raise ValueError('unclear the type of name required')

    def get_xs_app(self, xs): #!# change to get xs_app
        """
        Default approximation is in xs units
        """
        return 'e',self._approx_dic[xs]['app']

    def get_coef(self,xs=False):
        if xs!=False:
            return self._approx_dic[xs]['coef']
        else:
            return self._approx_dic

    def do_I_belong(self,taser_dict):
        """
        check if the key/vals of the dict corresponde to content of the class. Doesn't mean its the same
        """
        ret_flag=True
        for key,val in taser_dict.iteritems():
            try:
                if self.recived_app[key]!=val:
                    #print key,self.recived_app[key],val
                    ret_flag=False
            except KeyError:
                #print 'this is not in class key=', key
                ret_flag=False
        #print 'here',taser_dict,ret_flag
        return ret_flag

    def set_result(self, kind, sub_kind, result, **kwargs):
        """
        Unified set result frame for xs, mxs evaluation and error
        I choose kind/unit/xs because xs gets consumed later by the df
        """
        if not isinstance(result, dict):
            raise ValueError('wrong eval type')

        dict_type = kwargs.pop('type')
        if kwargs:
            raise ValueError('wrong kwargs')

        if dict_type == 'eval':
            dict_h = self._evaluation_dict
        if dict_type == 'obs':
            dict_h = self._obs_dict
        if dict_type == 'error':
            dict_h = self._error_dict

        for unit, data in result.iteritems():
            dict_h = self.initialize_dict(dict_h, kind, unit)
            if sub_kind in dict_h[kind][unit]:
                raise RuntimeError('Over-writing data in %s,%s' % (unit, sub_kind))
            dict_h[kind][unit][sub_kind] = data

    @staticmethod
    def initialize_dict(dic, kind, unit):
        """
        initializes dict
        """
        if kind not in dic:
            dic[kind] = OrderedDict()
        if unit not in dic[kind]:
            dic[kind][unit] = OrderedDict()
        return dic

    def dict_2_multindex_df(self, data_dict, nametuple_space, names, df_y_name, unit_flag=None):
        """
        this works for xs or error dict and for xs,Mxs,k...
        the problem is that IDK how to plot it properly
        """
        data_lst = []
        index_lst = []
        for key, tup in nametuple_space.iteritems():
            if unit_flag is None:
                data_lst.append(data_dict[key])  # for eval dict
            else:
                data_lst.append(data_dict[key][unit_flag])  # for error dict
            index_lst.append(list(tup))

        column_lst = list(np.transpose(np.array(index_lst)))
        data = np.transpose(np.array(data_lst))
        df = pd.DataFrame(data, columns=column_lst)
        df.columns.names = names
        df.index.name = df_y_name

        return df

    def dict_2_flat_df(self, data_dict, nametuple_space, IS_names, unit, df_y_name):
        """
        Tranforming to a FLAT dataframe. this is easer to hande in seaborn.
        Consider a vertical attachement of the different XS and thus repeat the grid-points many times. I'm not shure that an horizontal attachment would work well
        with seaborn i.e. having many columns with the same name.. In any case the ideal would be a multiindex with no repetition. But for that the seabonrn features must be somehow executed on the indexing of the dataframe and IDK how to od that.
        """

        vals = zip(*list(self.get_prp('test','grid_nmlz_val')))

        aux_dict = defaultdict(list)
        aux_grid = zip(zip(*self.get_prp('test','grid_nmlz_val')),
                       self._supp_inst.get_prp_valid('input_space_names',v_bypass=True))
        for key, tup in nametuple_space.iteritems():
            data = data_dict[unit][key]
            aux_dict[df_y_name].append(data)
            # this gets repeated len(OS) times in a column-flatened-df
            for opened_grid, IS_name in aux_grid:
                aux_dict[IS_name].append(list(opened_grid))

            # adding repeted names of OS
            for OS_name, OS_var in zip(tup._fields, tup):
                aux_dict[OS_name].append([OS_var for _ in xrange(len(data))])

        for key, val in aux_dict.iteritems():
            # print key, len(val), len(val[0]), len(list(itertools.chain(*val)))
            aux_dict[key] = list(itertools.chain(*val))

        # it would be nice to define the order better
        df = pd.DataFrame.from_dict(aux_dict, orient='columns')
        return df

    def get_err_data_multiindex(self, **kwargs):
        """
        Not used currently but it should, as soon as the interaction Multiindex/seaborn is solved
        """

        kind = kwargs['kind']
        sub_kind = kwargs['sub_kind']
        unit = kwargs['unit']

        # print self.OS_kind_names[kind]
        s_h = slice(None)
        template = [s_h for _ in self.OS_kind_names[kind]]
        if sub_kind == 'all':
            data_names = [slice(None)]  # nothing is sliced, duh
            #!# se supone que tiene q ser 'g'
            index_kind_output = self.multi_index_order[kind]['g']
        if sub_kind == 'g-wise':
            data_names = self._xirg['G']['name']
            index_kind_output = self.multi_index_order[kind]['g']
        if sub_kind == 'r-wise':
            data_names = self._xirg['R']['name']
            index_kind_output = self.multi_index_order[kind]['r']
        if sub_kind == 'i-wise':
            data_names = self._xirg['I']['name']
            index_kind_output = self.multi_index_order[kind]['i']

        df_vecs = []
        for xirg in data_names:
            aux = template
            aux[index_kind_output] = xirg
            df_vecs.append(tuple(aux))
        if not df_vecs:
            raise ValueError('Unable to select from the DataFrame')

        if any([not isinstance(xirg, str) for xirg in data_names]):  # for catching 'all basically'
            if all([isinstance(xirg, slice) for xirg in data_names]) and len(data_names) == 1:
                data_names = ['all']
            else:
                raise ValueError('dont know the generated sub_kind space')
        return data_names, [self._error_df_holder[kind][unit].loc[:, df_vec] for df_vec in df_vecs]

    def get_err_data_flat(self, **kwargs):
        """
        Obtaining sliced df for a statistical mission

        in general integrated over all IS and OS except that specifically indicated
        """
        consider_MACRT = None
        kind = kwargs.pop('kind')
        sub_kind = kwargs.pop('sub_kind')
        unit = kwargs.pop('unit')
        try:
            consider_MACRT = kwargs.pop('consider_MACRT')
        except KeyError:
            pass

        if kwargs:
            raise ValueError('wrong kwargs', kwargs)

        def extract(df):
            return df.loc[:, 'error']

        df = self._error_df_holder[kind][unit]

        if consider_MACRT == 'no' and kind == 'xs':
            df = df[~df['i'].isin(['MACRT'])]  # eliminating MACRT from DF. Default behaivour

        if consider_MACRT:
            if kind == 'mxs':
                raise ValueError('this shouldnt be here')

        if sub_kind == 'all':
            return [sub_kind], [extract(df)]
        if sub_kind == 'g-wise':
            data_names = self._xirg['G']['name']
            index_kind_output = 'g'
        if sub_kind == 'r-wise':
            data_names = self._xirg['R']['name']
            index_kind_output = 'r'
        if sub_kind == 'i-wise':
            data_names = self._xirg['I']['name']
            index_kind_output = 'i'
            data_names = self._xirg['I']['name']
        if sub_kind == 'xs-wise':
            """
            all the processing that follows could be unified with the pressiding in a single framework
            this clearly surpasess the other subkinds and they can be included here
            """
            get_OS = []
            order_aux = []
            # !# cant just use self.OS_kind_names[kind] ?
            for os_i in [cl for cl in df.columns if cl in self.OS_kind_names[kind]]:
                # get_OS.append({os_i: getattr(df, os_i).unique()})
                get_OS.append(getattr(df, os_i).unique())
                order_aux.append(os_i)
            sub_kinds = []
            for el in list(itertools.product(*get_OS)):
                sub_kinds.append({key: val for key, val in zip(order_aux, el)})
            ret_df = []
            data_names = []
            for sub_kind in sub_kinds:
                dfaux = df
                for col, val in sub_kind.iteritems():
                    dfaux = dfaux[dfaux[col] == val]
                dfaux = extract(dfaux)
                if dfaux.shape[0] != 0:
                    ret_df.append(dfaux)
                    data_names.append('_'.join([sub_kind[var] for var in self.OS_kind_names[kind]]))
            return data_names, ret_df

        return data_names, [extract(df[(df[index_kind_output] == var)]) for var in data_names]

    def save_sts_error(self, data, **kwargs):
        """
        Saves statistical results in dict. This should be eliminated in a future for working with df only
        """

        kind = kwargs.pop('kind')
        sub_kind = kwargs.pop('sub_kind')  # g-wise
        data_name = kwargs.pop('current_data')  # ['1','2']
        unit = kwargs.pop('unit')
        consider_MACRT = kwargs.pop('consider_MACRT', None)
        if kwargs:
            raise ValueError('bad kwargs')

        # this is too imperative. it should be only one flat loop on path
        if kind not in self._statistics_dict:
            self._statistics_dict[kind] = OrderedDict()

        if sub_kind not in self._statistics_dict[kind]:
            self._statistics_dict[kind][sub_kind] = OrderedDict()

        if data_name not in self._statistics_dict[kind][sub_kind]:
            self._statistics_dict[kind][sub_kind][data_name] = OrderedDict()

        if unit not in self._statistics_dict[kind][sub_kind][data_name]:
            self._statistics_dict[kind][sub_kind][data_name][unit] = OrderedDict()

        if not isinstance(data.values()[0], (list, float)):
            raise ValueError('Nested dict detected, only flat data here')

        for sts, val in data.iteritems():
            self._statistics_dict[kind][sub_kind][data_name][unit][sts] = val

        # statistical conversion
    def consolidate_stc_df(self):
        aux_dict = {}
        for key,val in dict(self.dict_common_app.items()+self.dict_special_app.items()).iteritems():
            aux_dict[key]=val

        _statistics_df_holder = OrderedDict()
        for kind in self._statistics_dict:
            _statistics_df_holder[kind] = OrderedDict()
            for sub_kind in self._statistics_dict[kind]:
                _statistics_df_holder[kind][sub_kind] = OrderedDict()
                for sub_kind_output in self._statistics_dict[kind][sub_kind]:
                    _statistics_df_holder[kind][sub_kind][sub_kind_output] = OrderedDict()
                    for unit in self._statistics_dict[kind][sub_kind][sub_kind_output]:
                        data = self._statistics_dict[kind][sub_kind][sub_kind_output][unit]
                        df = pd.DataFrame()
                        df = pd.DataFrame(data, index=[0])
                        df = df.assign(**aux_dict)
                        _statistics_df_holder[kind][sub_kind][sub_kind_output][unit] = df

        self._statistics_df_holder = _statistics_df_holder
        # this could be done automaticly using some path variable
        self._reconstruct_order_sttc = ['kind', 'sub_kind']

    def consolidate_data_df(self, kind_vec=None):
        """Transofrmation of dicts into pd. I hope to end up working in pd only.
        #!# I don't like very much the dinamic creation getattr(self, kind + '_nametuple_space') as the exposed property appears everywhere and is kind of brittle.
        #!# this should be unified with the sts consolidation, though that may not be tirivial...
        """
        # producing a FG for the multiindex space. Multiindex is kind of complicated
        if kind_vec is None:
            kind_vec = self._present_kinds
        transformation = self.dict_2_flat_df

        self._error_df_holder = OrderedDict()
        self._eval_df_holder = OrderedDict()
        self.multi_index_order = OrderedDict()
        for kind in kind_vec:
            names = self.OS_kind_names[kind]
            self.multi_index_order[kind] = {key: pos for key,
                                            pos in zip(names, range(len(names)))}
            self._error_df_holder[kind] = OrderedDict()
            self._eval_df_holder[kind] = OrderedDict()
            for unit in self._error_dict[kind].keys():
                self._error_df_holder[kind][unit] = transformation(
                    self._error_dict[kind], self.get_prp('a2', kind + '_nametuple_space'), names, unit, 'error')

            unit = self._evaluation_dict[kind].keys()
            if len(unit)>1:
                raise RuntimeError('there shouldnt be more than one unit')
            unit=unit[0]
            self._eval_df_holder[kind][unit] = transformation(
                self._evaluation_dict[kind],self.get_prp('a2',kind + '_nametuple_space'), names, unit, 'eval')

        self._reconstruct_order_eval = ['kind', 'unit']
        self._reconstruct_order_error = ['kind', 'unit']

    def get_df_holder(self, typE, **path):
        """
        return the holder of the df, only the first two level of df are required in path
        example: app.get_df_holder('error', kind=type_df, unit='e')
        A more comprehensive ando robust system could be developed later, for geting all dfs
        """
        h = None
        if typE == 'eval':
            h = self._eval_df_holder
            reco = self._reconstruct_order_eval
        if typE == 'error':
            h = self._error_df_holder
            reco = self._reconstruct_order_error
        if typE == 'sttc':
            h = self._statistics_df_holder
            reco = self._reconstruct_order_sttc
        if not h:
            raise ValueError('could not find the df ', typE)
        if path:
            reconstructed_path = [path.pop(key) for key in reco if key in path]
            if path:
                raise ValueError('bad kwargs')
            return self.nested_access(h, *reconstructed_path)
        return h

    def nested_access(self, data, *args):
        """
        Allows for reciving a dic and returns elements in the path defined ar args.
        If there is non, None is returned
        This facilitates  conditional path definitions as allows helps reducing notation for long paths
        Usage:
        dic_path=['A2','UO2 Dm=0.72 Tf=600','data','I','U235']
        print nested_access(apollo2_data, *dic_path).keys()
        """
        if args and data:
            element = args[0]
            if element:
                value = data.get(element)
                return value if len(args) == 1 else self.nested_access(value, *args[1:])

    def del_dicts(self):
        """
        Deleting auxiliary dict to save only df. The difference in final size is big
        """
        del(self._statistics_dict)
        del(self._error_dict)
        del(self._evaluation_dict)

    def approximate_xs(self,user_xs_names=False):  # this parent?
        """
        xs approximation is built for every xs. A dict of app and properties is returned, this may change for different types of approximations
        """
        _approx_dic = OrderedDict()
        if user_xs_names:
            xs_names=user_xs_names
        else:
            xs_names=self.get_prp('a2','xs_namespace')
        for xs_name in xs_names:
            _approx_dic[xs_name] = {}
            xs = list(self.get_fprp('supp','xs',kind=xs_name))
            ret_dict=self.represente(xs,self._supp_inst.cartesian_val,self.get_prp('supp','grid_nmlz_val',sub_kind=xs_name),xs_name=xs_name)
            for ret_n,ret_v in ret_dict.iteritems():
                _approx_dic[xs_name][ret_n]=ret_v

        self._approx_dic = _approx_dic

class Approximation_SPL(Approximation_Container):
    """Spline interpolation class

    Responsabilities:
    *build approximation(coeficients)
    *evaluate
    """
    __doc__ = API()

    # go to parent
    def __init__(self,app_n, app, supp):
        super(Approximation_SPL, self).__init__(supp,app_n,app['supp_type'])
        self.recived_app=app
        self._orders = app['orders']
        try:
            self._degree = [val-1 for val in app['orders']] 
        except TypeError:
            self._degree=['dummy']
        self._types = app['types']
        self._t = app['t']
        self.dict_special_app={}
        self.dict_special_app['degree']=self._degree[0]
        self.dict_special_app['t']=app['t'][0]
        self.dict_special_app['app_name']='f'+app['types'][0]+app['t'][0]+'dg'+str(self._degree[0])+'_'+app['supp_type']+'n'+str(supp.size)
        self.dict_special_app['app_name']=self.dict_special_app['app_name'].replace('__','_')
        self.dict_special_app['mth_name']='f'+app['types'][0]+app['t'][0]+'dg'+str(self._degree[0])+'_'+app['supp_type']
        self.dict_special_app['mth_name']=self.dict_special_app['mth_name'].replace('__','_')

    def represente(self, xs,tau_gen,supp_grid,xs_name=None):
        """Produce spline representation"""
        k = self._orders
        d = len(self._orders)
        type_aux = [s.replace('fg_bsp', '').replace('_', '') for s in self._t]
        if self._types[0]=='_bsp_':
            t = knot_vector(tau_gen, k, type_aux)  # define knot vector
            base_i = [S(k[i], t[i]) for i in range(d)]  # define based
        elif self._types[0]=='_pp_':
            base_i = [PP(len(tau_gen[i]), tau_gen[i]) for i in range(d)]  # define based
            self._orders=[base_i[i].k for i in range(d)]
            self._degree = [val-1 for val in self._orders]  # !# not clear
        else:
            raise ValueError('dont know how to repesent')
        if base_i is None:
            raise ValueError('something whent wrong when computing the base')
        Fd1 = Fd(base_i[0])
        for i in range(d)[1:]:
            Fd1 = Fd1 * base_i[i]
        fd1 = fd(Fd1)
        try:
            fd1.cmpcoef(tau_gen, xs)
        except ValueError:
            raise ValueError
        ret_dic={'app':fd1}
        if self._types[0]=='_bsp_':
            ret_dic['t']=t
        else:
            ret_dic['t']=tau_gen
        try:
            ret_dic['coef']=np.concatenate(fd1.coef)
        except ValueError:
            ret_dic['coef']=fd1.coef
        return ret_dic

    def __str__(self):
        return 'spline' + str(self._types) + str(self._orders)

class Approximation_CKM(Approximation_Container):
    """
    Kernel method approximation using spline kernel developed in fortran

    to build module from C
    swig -python kernels.i
    python setup.py build_ext --inplace
    Apprioxmiation C Kernel Method

    FORTRAN
    python2.7 -m numpy.f2py kernels.f90 -m fortran_k -h fortran_k2.pyf --overwrite-signature
    python2.7 -m numpy.f2py --opt=-O2  -c fortran_k2.pyf kernels.f90

    cProfile
    python -m cProfile -o AL xs_main_building_app.py
    python -m pstats AL
    help, sort, sort cumulative,...
    tottime very usefull
    stats 20


    """

    def __init__(self,app_n, app, supp):
        super(Approximation_CKM, self).__init__(supp,app_n,app['supp_type'])
        # how to take different type and amount of specification for the kernel could be managed in a better way
        self.recived_app=app
        self.app_name='f'
        self._types = app['types']
        self.app_name=self.app_name+app['types'][0]
        self._k = app['k'][0]
        self.app_name=self.app_name+app['k'][0]
        try:
            self._degree = app['degree']
            self._orders = [a+1 for a in app['degree']]
            self.app_name=self.app_name+'_dg'+str(app['degree'][0])+'_'
        except KeyError:
            pass
        try:
            self._omega = app['omega']
            self.app_name=self.app_name+'_om'+str(app['omega'][0])+'_'
        except KeyError:
            pass
        try:
            if app['norm']=='_ln(bu)_':
                self._ISnorm_v=np.log
                self._ISnorm=app['norm']
                self.app_name=self.app_name+self._ISnorm
            else:
                raise ValueError('cant find norm')
        except KeyError:
            self._ISnorm=''
            self._ISnorm_v=None
        self.app_name=self.app_name+app['supp_type']
        self.check_mode_flag=False
        try:
            self._regul=app['regul'][0]
        except KeyError:
            self._regul=False
        try:
            self._rcond=app['rcond'][0]
            self.app_name=self.app_name+'_rc'+str('{0:1.0e}'.format(self._rcond)).replace('1e-','')+'_'
        except KeyError:
            self._rcond=False
        try:
            self._precd=app['precd']
            self.app_name=self.app_name+self._precd
        except KeyError:
            self._precd=False
        try:
            self._reduceC=app['coef_reduced']
            self.app_name=self.app_name+'_RC'+str(self._reduceC)+'_'
        except KeyError:
            self._reduceC=False
        try:
            self._ALTOL=app['TOL']
            self.app_name=self.app_name+'_ALT'+str(self._ALTOL)+'_'
        except KeyError:
            self._ALTOL=False
        #confiando en el nombre..
        self.multiple_xs_valid_iterators=False
        self._K={}
        self._invK={}
        self.supp_conditioned_X={}
        self.len_X={}
        if 'wal_' in app['supp_type']:
            self.multiple_xs_valid_iterators=True
        else:
            pass
        self.calculation_of_invK=False

        def build_dict_special_app1(k):
            dic={}
            dic['types']=app['types'][0]
            dic['k']=app['k'][0]
            dic['app_name']=self.app_name+'n'+str(supp.size)+'_'
            dic['app_name']=dic['app_name'].replace('__','_')
            dic['mth_name']=self.app_name
            dic['mth_name']=dic['mth_name'].replace('__','_')
            return dic

        self.kernel=None
        if app['k'][0]=='_bn_':
            self._parameter1=self._degree
            self.kernel=KM._KernelBernoulli
            self.F_flag='KBN'

        if app['k'][0]=='_sp_':
            self._parameter1=self._degree
            self.kernel=KM._KernelSpline
            self.F_flag='KSP'

        if app['k'][0]=='_pb_':
            self._parameter1=self._omega
            self.kernel=KM._KernelPB
            self.F_flag='Not implemeted'

        self.dict_special_app=build_dict_special_app1(app['k'][0])
        if not self.kernel:
            print self.kernel
            raise ValueError('cant define kernel')
        self.acceleration_n=False


    def accelerate(self,acceleration):
        # like the test values to calculate k(x,xi) only once
        if 'pre_calc_K(x,xi)' in acceleration:
            self.acceleration_n='pre_calc_K(x,xi)'
            self.conditioned_X_test=self.IS_normalization(np.array(acceleration['pre_calc_K(x,xi)'],order='F'),self._precd,x_type='vec')
            self.K_test={}
        else:
            raise ValueError('acceleration invoked without understandable values')
        if self.F_flag==False:
            raise ValueError('acceleration presuposes use of Frotran functions')
        if self.multiple_xs_valid_iterators==True:
            raise ValueError('accelration for varying coef not implemeted')

        if self.acceleration_n=='pre_calc_K(x,xi)':
            xs_n=self.supp_conditioned_X.keys()
            if len(xs_n)>1:
                raise ValueError('acceleration does not admit changing alpha')
            xs_n=xs_n[0]
            def K_test_calcul_shell(xs_n):
                return FKM.calcul_k_generic(self.F_flag,self.len_flag,self._parameter1,self.supp_conditioned_X[xs_n],self.conditioned_X_test,self.len_k,self.len_X[xs_n],len(self.conditioned_X_test))
            for xs_n in self.supp_conditioned_X.keys():
                self.K_test[xs_n]=K_test_calcul_shell(xs_n)

    def represente(self, xs,tau_gen,supp_grid,xs_name=None):
        r"""Kernel method
        Let f=sum_{i=1}^{i=N}\alpha_i*k(x,x_i)
        """
        if xs_name not in self._K:
            if len(self._K.keys())==0 or self.multiple_xs_valid_iterators==True:
                try:
                    X=self.supp_conditioned_X[xs_name]
                except KeyError:
                    X=np.array(list(supp_grid),order='F')
                    X=self.IS_normalization(X,self._precd,x_type='vec')
                    self.supp_conditioned_X[xs_name]=X
                    self.len_X[xs_name]=int(len(X))
                try:
                    self.len_flag
                except AttributeError:
                    self.len_flag,self.len_k,=len(self.F_flag),int(len(self._parameter1))
                K,invK=self.calcul_invK(X,self._parameter1,len(self._parameter1),self.K_prod,self.kernel,self.len_flag,self.len_k,self.len_X[xs_name],F_flag=self.F_flag,regul=self._regul,rcond=self._rcond)
                self._K[xs_name]=K
                self._invK[xs_name]=invK
            else:
                invK=self._invK[self._K.keys()[0]]

        xsmax,norm_xsmean,xs_norm=self.normalize1(xs,safe=False)
        ret_dict={}
        alpha=np.array(np.matmul(invK,xs_norm),order='F')
        if self._reduceC!=False:
            alpha=self.coefficient_reduction(alpha,self._reduceC)
        ret_dict['coef']=alpha
        if self.check_mode_flag==True:
            try:
                self.check_norm_dict[xs_name]=xs_norm
            except AttributeError:
                self.check_norm_dict=OrderedDict()
                self.check_norm_dict[xs_name]=xs_norm
        ret_dict['xsmax']=xsmax
        ret_dict['norm_xsmean']=norm_xsmean
        return ret_dict

    def get_xs_app(self, xs):
        """
        the approximation 'app' is built on the fly but never stored since I had allot of problems with pickeling it
        The lighter the f is the better, here is where all filthy stuff must take place, and not inside f!
        I use val normailzation to avoid problems with the map(g,x)
        """

        X=self._aux_multi_xs(self.supp_conditioned_X,key=xs)
        len_X=self._aux_multi_xs(self.len_X,key=xs)
        norm_xsmean, xsmax=self._approx_dic[xs]['norm_xsmean'],self._approx_dic[xs]['xsmax']
        if self.F_flag:
            if self.acceleration_n=='pre_calc_K(x,xi)':
                """
                this dual behaivour allows for not modifing the analysis part
                """
                K_test=self._aux_multi_xs(self.K_test,key=xs)
                def f(x,idx=None):
                    """
                    if idx is given, as in the supervised learning then the approximation can be solved as an inner product of alpha and the precalculated K_test matrix
                    Otherwise an independend k(x,xi) calculation is performed
                    """
                    if idx!=None:
                        return self.Finner(idx,self._approx_dic[xs]['coef'],norm_xsmean,xsmax,K_test)
                    else:
                        x=self.IS_normalization(x,self._precd,x_type='point')
                        return self.Fpredict(self._parameter1,x,X,self._approx_dic[xs]['coef'],norm_xsmean,\
             xsmax,self.len_flag,self.len_k,len_X,self.F_flag)
            else:
                def f(x,idx=None):
                    x=self.IS_normalization(x,self._precd,x_type='point')
                    return self.Fpredict(self._parameter1,x,X,self._approx_dic[xs]['coef'],norm_xsmean,\
             xsmax,self.len_flag,self.len_k,len_X,self.F_flag)
        else:
            # no longer maintained
            def f(x,idx=None):
                x=self.IS_normalization(x,self._precd,x_type='point')
                return self.pypredict(self._parameter1,x,X,self.kernel,self._approx_dic[xs]['coef'],norm_xsmean,\
             xsmax,self.K_prod)
        return 'e',f

    def _aux_multi_xs(self,Dict,key=None):
        if self.multiple_xs_valid_iterators==False:
            return Dict[Dict.keys()[0]]
        else:
            return Dict[key]

    def IS_normalization(self,x,key,x_type=None):
        """
        has a vector as input
        """
        try:
            bu_idx_py=self.bu_idx_py
        except AttributeError:
            bu_idx_py=[tup.py_idx for tup in self.get_prp('a2','tuple_index_order') if tup.dim_name=='BURNUP'][0] #will be 0 in general
            self.bu_idx_py=bu_idx_py

        g=self.normalize2(key,bu_idx_py=bu_idx_py)
        x=list(cp.deepcopy(x))
        if x_type=='point':
            working_array=g(x) # for f evaluations
        if x_type=='vec':
            working_array=map(g,x) # for support preconditioning in coefficent obtantion
        if x_type==None:
            raise ValueError('type of x unknown')

        if key=='_sqrt(bu)_':
            if x_type=='vec':
                for xi,sqrtbui in zip(x,working_array):
                    xi[bu_idx_py]=sqrtbui
                    working_array=x
            if x_type=='point':
                x[bu_idx_py]=working_array
                working_array=x
        return np.array(working_array,order='F')

    @staticmethod
    def normalize2(key,bu_idx_py=None):
        """
        if no normalization the the identity is returned
        """
        if key=='_sqrt(all)_':
            return np.sqrt

        if key=='_sqrt(bu)_':
            return lambda xi: np.sqrt(xi[bu_idx_py])

        return lambda xi:xi

    @staticmethod
    def normalize1(xs, safe=True):
        """
        normalization of OS
        Required normalization for approximating with KM
        If allot of time is detected in this step normalized XS can be stored
        """
        if safe==True:
            if min(xs)<0:
                raise ValueError('normilizing negative XS')
            if type(xs)!=list:
                raise ValueError('this was usually a list')
        xsmax=1
        gtau=np.array([g/xsmax for g in xs],order="C")
        norm_xsmean=np.mean(gtau)
        return xsmax,norm_xsmean,np.array([g-norm_xsmean for g in gtau],order="C")


    def coefficient_reduction(self,alpha,percentage):
        """
        instead of actually eliminating the coefficient I just put it at 0 as no not need to modify enything else in the class.
        Then the size of the supprot is just the non-zero coefficients
        """
        if percentage<0 or percentage>1:
            raise ValueError('limits of reduction dont make sense')

        indexed_alpha=list(enumerate(abs(alpha)))
        ordered_indexed_alpha=list(reversed(sorted(indexed_alpha,key=lambda x:x[1])))
        retianed_index=[val[0] for val in ordered_indexed_alpha[0:int(len(ordered_indexed_alpha)*(1-percentage))]]
        if len(retianed_index)<1:
            raise ValueError('not enough coefficients could be retained')
        ret_alpha=cp.deepcopy(alpha)
        for idx,_ in enumerate(ret_alpha):
            if idx not in retianed_index:
                ret_alpha[idx]=0.0
        self.dict_common_app['size']=np.count_nonzero(ret_alpha) #note that this is performed for every XS
        return ret_alpha
        

    @staticmethod
    def calcul_invK(X,k,d,Kprod,kernel,len_flag,len_k,len_X,F_flag=False,regul=False,rcond=False):
        def pyK():
            K=np.zeros([len(X),len(X)],order="C")
            for i,x in enumerate(X):
                for j,y in enumerate(X):
                    K[i][j]=Kprod(x,y,k,kernel)
            return K
        def FK():
            return FKM.calcul_k(F_flag,len_flag,X,k,len_k,len_X)

        if F_flag!=False:
            K=FK()
        else:
            print 'warning using python insted of F'
            K=pyK()

        if regul!=False:
            K=K +regul*len(K)*np.eye(len(K))
            print regul*len(K)*np.eye(len(K))
            print regul
            invK=np.linalg.inv(K)
            sys.exit()
        else:
            if rcond==False:
                rcond=1E-15
            invK=np.linalg.pinv(K,rcond=rcond)
        return K,invK


    @staticmethod
    def Finner(idx,alpha,norm_xsmean, xsmax,K_test):
        """
        Attention coherence of searchin the right index must be assured by construction. i.e. in the generation of K_test the busequent order of idx must not change
        """
        ret_val=FKM.inner_product(alpha,K_test[:,idx],len(alpha))
        return (ret_val+norm_xsmean)*xsmax

    @staticmethod
    def pypredict(k,x,X,kernel,alpha,norm_xsmean, xsmax,Kprod):
        """
        Making the predcition a static method allows to use it outside the class
        """
        ret_val=np.inner(alpha,[Kprod(x,xi,k,kernel) for xi in X])
        return (ret_val+norm_xsmean)*xsmax

    @staticmethod
    def Fpredict(k,x,X,alpha,norm_xsmean, xsmax,len_flag,len_k,len_X,F_flag):
        ret_val=FKM.predict(F_flag,len_flag,X,alpha,x,k,len_k,len_X)
        return (ret_val+norm_xsmean)*xsmax


    @staticmethod
    def K_prod(x,y,k,kernel,p=1):
        for xi,yi,ki in zip(x,y,k):
            p=p*kernel(xi,yi,ki)
        return p


    def py_Bernoulli(self,x,degre):
        """
        Trying to solve BN1 problem
        """
        y=None
        if degre == 0: y = 1.
        if degre == 1: y = x - 0.5
        if degre == 2: y = x * x - x + 1./6.
        if degre == 3: y = x * (x * (x - 1.5) + 0.5)
        if degre == 4: y = x * (x * (x * (x - 2.0) + 1.0)) - 1./30.
        if degre == 5: y = x * (x * (x * (x * (x - 2.5) + 5./3.)) - 1./6.)
        if degre == 6: y = x * (x * (x * (x * (x * (x - 3.0) + 2.5)) - 0.5)) + 1./42.
        return y

    def py_KernelBernoulli(self,x, y, degre):
        #double v, s, tgamma(double), pow(double, double), floor(double);
        s=0
        for d in range (degre+1):
            s = s + self.py_Bernoulli(x, d) * self.py_Bernoulli(y, d) / mth.factorial(d)**2
        v=abs(x-y)# parte entera o abs?
        #v=x-y-mth.floor(x-y)
        return s + ( ((-1)**(degre+1.)) /mth.factorial(2. * degre ))* self.py_Bernoulli(v, 2 * degre)


    def __str__(self):
        return str( self.dict_special_app['app_name'])


class Approximation_SK(Approximation_Container):
    """
    Kernel method approximation using scikit learn gaussian kernels
    """
    def __init__(self,app_n,app, supp):
        super(Approximation_SK, self).__init__(supp,app_n,app['supp_type'])
        self._types = app['types']
        self.dict_special_app={}
        self.dict_special_app['types']=app['types'][0]
        self.dict_special_app['app_name']='f'+app['types'][0]+app['k'][0]+app['supp_type']+'_'+'n'+str(supp.size)
        self.dict_special_app['app_name']=self.dict_special_app['app_name'].replace('__','_')
        self.dict_special_app['mth_name']='f'+app['types'][0]+app['k'][0]+app['supp_type']
        self.dict_special_app['mth_name']=self.dict_special_app['mth_name'].replace('__','_')

        param_grid={}
        param_grid["alpha"]= [1e-2,1e-3,1e-4,0]
        # param_grid["alpha"]= [0]
        if '_sg_' in app['k']:
            param_grid['kernel']=['sigmoid']
            param_grid['gamma']=[0.1,10,100]
            param_grid['coef0']=[0.1,1,10,100]

        if '_poly_' in app['k']:
            param_grid['kernel']=['polynomial']
            param_grid['gamma']=[0.1,10,100]
            param_grid['coef0']=[0.1,1,10,100]
            param_grid['degree']=[2,3,4]

        if '_lap_' in app['k']:
            param_grid["kernel"]=['laplacian']
            param_grid['gamma']=[0.1,10,100]

        if '_ln_' in app['k']:
            param_grid["kernel"]=['linear']

        if '_rbf_' in app['k']:
            param_grid["kernel"]=['rbf']
            param_grid['gamma']=[0.1,10,100]

        if '_mt_' in app['k']:
            from sklearn.gaussian_process.kernels import Matern
            a_vec=[0.1, 1,10]
            b_vec=[0.5,1.5,2.5]
            param_grid["kernel"]=[Matern(length_scale=a,nu=b) for a in a_vec for b in b_vec]

        if '_ch_' in app['k']:
            param_grid["kernel"]=['chi2']
            param_grid['gamma']=[0.1,10,100]

        self._km=GridSearchCV(KernelRidge(), cv=5,param_grid=param_grid,refit=True,scoring='r2',verbose=0) #,refit=True,scoring='r2'


    @staticmethod
    def normalize1(xs):
        """
        Required normalization for approximating with KM
        If allot of time is detected in this step normalized XS can be stored
        """
        if min(xs)<0:
            raise ValueError('normilizing negative XS')
        xsmax=max(xs)
        gtau=np.array([g/xsmax for g in xs],order="C")
        norm_xsmean=(np.max(gtau)+np.min(gtau))/2
        return xsmax,norm_xsmean,np.array([g-norm_xsmean for g in gtau],order="C")

    def represente(self, xs,tau_gen,supp_grid,xs_name=None):
        """Kernel method
        Let f=sum_{i=1}^{i=N}\alpha_i*k(x,x_i)
        """
        X=list(supp_grid)

        xsmax,norm_xsmean,xs_norm=self.normalize1(xs)
        ks=cp.deepcopy(self._km)
        ks.fit(np.array(X),xs_norm)
        ret_dict={}
        ret_dict['xsmax']=xsmax
        ret_dict['norm_xsmean']=norm_xsmean
        ret_dict['fitted']=ks
        ret_dict['score']= ks.best_score_
        ret_dict['param']= ks.best_params_
        return ret_dict

    def get_xs_app(self, xs):
        """
        the approximation 'app' is built on the fly but never stored since I had allot of problems with pickeling it
        """
        def f(x):
            norm_xsmean, xsmax=self._approx_dic[xs]['norm_xsmean'],self._approx_dic[xs]['xsmax']
            ret_val=self._approx_dic[xs]['fitted'].predict(np.array(x).reshape(1, -1))[0] #it returns a list with one element

            return (ret_val+norm_xsmean)*xsmax
        return 'e',f

    def get_coef(self,xs):
        return 'to implement'

    def __str__(self):
        return 'scikit-learn' + str(self._types)  + str(self._supp_inst.size)

