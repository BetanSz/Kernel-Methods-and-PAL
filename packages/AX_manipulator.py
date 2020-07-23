import h5py
import subprocess
from packages.tools_SG_common import *
from packages.tool_FP import *

"""
Procedural aproach of data transofrmation.
Al a lower level an imperative apporach is necesarry, as procedures must be customized to incomming
 data. However FP paradigm has been used to sme extent. Higher order fucntions are presents,
 functions are pure, and a procedure in general overwrigths original data, emulating the type
 of data transofmration seen with inmutable data.
"""


class XsTupPickable(namedtuple('xs_tup', 'x i r g')):
    pass


class StateParPickable(namedtuple('grid_point', 'dim_name py_idx a2_idx')):
    pass


def standarize_apollo2data(apollo2_data, num_dim, retained_i, retained_r, retained_g, verbose=False,res_nufi_removal=True):
    """FP-style of data treatment for a2 for the function space. Standarization in-situ, i.e.
     without coping the object"""

    if verbose == True:
        print ' ... standarization of IS and OS'
    apollo2_data = standarization_oftags(apollo2_data)
    apollo2_data = standarization_ofxs(
        apollo2_data,res_nufi_removal, retained_i=retained_i, retained_r=retained_r, retained_g=retained_g)
    # FP-style for the phase-space
    apollo2_data['PS'] = standarization_ofPS(
        apollo2_data['PS'], num_dim, apollo2_data['order_tuple'])

    try:
        apollo2_data.pop('order_tuple')
    except KeyError:
        warnings.warn('Check structure of incomming a2 dic')

    if verbose == True:
        print ' ... standarization of conc'

    apollo2_data['PS']['real']['conc'] = standarization_ofconc(apollo2_data['I'])
    apollo2_data['I'] = elimination_ofconc(apollo2_data['I'])
    try:
        apollo2_data['fi'].pop('iota')
        apollo2_data['k'].pop('b2')
    except KeyError:
        warnings.warn('Check structure of incomming a2 data. dic.keys() = %s' % apollo2_data.keys())

    # FP-style for the k/SI
    if verbose == True:
        print ' ... standarization of k,fi'
    for eigen_pair in ['k', 'fi']:
        for eigen in apollo2_data[eigen_pair]:
            if eigen != 'b2':  # this is not always present
                if verbose == True:
                    print ' ... doing ', eigen_pair, eigen
                tmp = standarization_ofPScodomain(
                    apollo2_data[eigen_pair][eigen], apollo2_data['PS'])
                apollo2_data[eigen_pair][eigen] = OrderedDict()
                apollo2_data[eigen_pair][eigen][eigen] = tmp
                if eigen_pair == 'k':
                    apollo2_data[eigen_pair][eigen]['ro'] = reactivity_calculation(
                        apollo2_data[eigen_pair][eigen][eigen])
    # FP-style for the xs data
    # easier to use auxiliary dic than in-situ transofrmation, but less elegant
    # handeling the interface at this level, thus pure functions can be used
    if verbose == True:
        print ' ... standarization of xs'
    aux_dic = OrderedDict()
    output_space_names = []
    output_space_tup = OrderedDict()
    for x in ['1']:  # eventually zones will come in a2 data
        aux_dic[x] = OrderedDict()
        for i in apollo2_data['I']:
            aux_dic[x][i] = OrderedDict()
            for r in apollo2_data['I'][i]['R']:
                aux_dic[x][i][r] = OrderedDict()
                for g in apollo2_data['I'][i]['R'][r]:
                    if verbose == True:
                        print ' ... doing ', i, r, g
                    aux_dic[x][i][r][g] = standarization_ofPScodomain(
                        apollo2_data['I'][i]['R'][r][g], apollo2_data['PS'])
                    name = '_' + x + '_' + i + '_' + r + '_' + g
                    output_space_names.append(name)
                    output_space_tup.update({name: XsTupPickable(x, i, r, g)})
    apollo2_data.pop('I')
    apollo2_data['xs'] = aux_dic
    apollo2_data['xs_names'] = output_space_names
    apollo2_data['xs_tuple'] = output_space_tup
    del(aux_dic)  # the reference is eliminated
    apollo2_data['XIRG'] = hierachy_nomenclature(apollo2_data)
    return apollo2_data


def hierachy_nomenclature(a2_data):
    """Returns the namespace and the cardinality of the output space"""
    ret_dic = OrderedDict()
    ret_dic['X'] = OrderedDict()
    ret_dic['X']['name'] = a2_data['xs'].keys()
    ret_dic['X']['N'] = len(a2_data['xs'].keys())
    ret_dic['I'] = OrderedDict()
    ret_dic['I']['name'] = a2_data['xs']['1'].keys()
    ret_dic['I']['N'] = len(a2_data['xs']['1'].keys())
    ret_dic['R'] = OrderedDict()
    ret_dic['R']['name'] = a2_data['xs']['1']['U235'].keys()
    ret_dic['R']['N'] = len(a2_data['xs']['1']['U235'].keys())
    ret_dic['G'] = OrderedDict()
    ret_dic['G']['name'] = a2_data['xs']['1']['U235']['abso'].keys()
    ret_dic['G']['N'] = len(a2_data['xs']['1']['U235']['abso'].keys())
    return ret_dic


def reactivity_calculation(a2_data):
    """Calculation of reactivity"""
    ret_val = np.array([1E5 * (point - 1) / point for point in a2_data])
    ret_val.setflags(write=False)
    return ret_val


def standarization_ofPScodomain(a2_data, PS):
    """
    for data comming in the form
    (tuple of instantaneus values):vector of values in increasing burnup
    the corresponding pppack form is found, corresponding the the grid indexation. This is done matching the order of the codomain
    withe the one on the grid.
    """

    if 'grid' not in PS['nmlz']['indx'].keys():
        raise RuntimeError('py-user indexed grid required for mapping codomain')

    # generating mapping from a2 to py in FP style
    # ATTENTION: "-1" is substracted from the return a2 tuple due to  a2 starting indexing with 1
    # while this info is not present in the incoming data this has to be hardcoded in the mapping
    map_a2_2_py = []
    for dimi in PS['real']['data'].keys():
        if dimi.dim_name.lower() == "burnup":
            if dimi.a2_idx is not None:
                raise RuntimeError('assumed indexation g(iota,bu) is false')
        else:
            map_a2_2_py.append(lambda x, i_bound=dimi.a2_idx: x[i_bound] - 1)

    # print map_a2_2_py
    # sys.exit()
    # flatting a2 data as x_{py-grid-point}, g(x)
    gx_unordered = []
    for a2_tuple in a2_data:
        if type(a2_tuple) != tuple:
            raise RuntimeError('cant understand codomain')
        for bu_idx, val in enumerate(a2_data[a2_tuple]):
            gx_unordered.append([tuple([bu_idx] + [f(a2_tuple) for f in map_a2_2_py]), val])

    grid_len = len(PS['nmlz']['indx']['grid'])
    if len(gx_unordered) != grid_len:
        print PS['nmlz']['indx']['grid']
        print gx_unordered
        raise RuntimeError('Mismatch between the length of a2 data and corresponding grid points')

    # ordering the calc points according to the py-grid (ppack order)
    gx_pygridindex_unordered = []
    for x, gx in gx_unordered:
        try:
            # stored as [required position, value]
            gx_pygridindex_unordered.append([PS['nmlz']['indx']['grid'].index(x), gx])
        except ValueError:
            print 'a2 point could not be found in grid'
            sys.exit()

    if len(gx_pygridindex_unordered) != grid_len:  # this should be Order(1)
        raise RuntimeError('Mismatch between the length of a2 data and corresponding grid points')
    gx_ordered = np.array([gx[1] for gx in sorted(
        gx_pygridindex_unordered, key=lambda x:x[0])], order='F')
    gx_ordered.setflags(write=False)
    return gx_ordered

def standarization_oftags(a2_data, tag0='UO2 Dm=0.72 Tf=600'):
    """
    Cheking only one JdD is present
    """
    if 'A2' in a2_data:
        tag = a2_data['A2'].keys()
        if len(tag) != 1:
            tag = tag0
            print 'warning, tags=', tag
        else:
            tag = tag[0]
        return a2_data['A2'][tag]['data']
    if 'data' in a2_data:
        return a2_data['data']
    raise RuntimeError('cant read a2 data')


def standarization_ofxs(a2_data,res_nufi_removal, **kwargs):
    """ eliminating un-wanted parts of the xs values, or general reduction """
    try:
        retained_r = kwargs.pop('retained_r')
    except KeyError:
        retained_r = a2_data['I'][a2_data['I'].keys()[0]]['R'].keys()
    try:
        retained_i = kwargs.pop('retained_i')
    except KeyError:
        retained_i = a2_data['I'].keys()
    try:
        retained_g = kwargs.pop('retained_g')
    except KeyError:
        retained_g = a2_data['I'][a2_data['I'].keys()[0]]['R']['abso'].keys()

    if kwargs:
        raise TypeError('Unexpected **kwargs: %r' % kwargs)
    a2_data = isotopes_remotion(a2_data, retained_i,residual_rebuild=False) # ir order to do residual_rebuild its necessary to include res nufi in calculations. Not happening for now..
    a2_data = reaction_remotion(a2_data, retained_r)
    a2_data = group_remotion(a2_data, retained_g)
    a2_data = nonphysicalxs_remotion(a2_data,res_nufi_removal)
    return a2_data


def standarization_ofPS(a2_data, num_dim, a2_tuple):
    """ Standarization of the PS"""
    py_order = ['BURNUP', 'TEMPERATURE_COMBUSTIBLE',
                'CONCENTRATION_EN_BORE', 'DENSITE_MODERATEUR', 'PISCINE',
                'BURNUPstep', 'CONTROL_ROD', 'POWER_DENSITY', 'PHASE']
    # interface is manege here to make functions as pure as possible
    # a2_data['PS'] = efective_real_state_space(a2_data['PS'], num_dim)
    temp = efective_real_state_space(a2_data, num_dim)
    a2_data = OrderedDict()
    a2_data['real'] = OrderedDict()
    a2_data['real']['data'] = ordered_real_state_space(temp, py_order, a2_tuple)
    a2_data['nmlz'] = OrderedDict()
    a2_data['nmlz']['tau_gen'] = normalize_state_space(a2_data['real']['data'])
    a2_data['nmlz']['N_tau'] = [len(vec) for vec in a2_data['nmlz']['tau_gen'].values()]
    a2_data['nmlz']['N_prod'] = np.prod(a2_data['nmlz']['N_tau'])
    a2_data['nmlz']['grid'] = to_pppack_style(a2_data['nmlz']['tau_gen'].values())
    if len(a2_data['nmlz']['grid']) != a2_data['nmlz']['N_prod']:
        raise RuntimeError('Generation of grid did not reproduce Cartesian construction %r != %r'
                           % (len(a2_data['nmlz']['grid']), a2_data['nmlz']['N_prod'])
                           )
    a2_data['nmlz']['indx'] = OrderedDict()
    a2_data['nmlz']['indx']['tau'] = [[idx for idx, _ in enumerate(vec)] for vec in a2_data[
        'nmlz']['tau_gen'].values()]
    a2_data['nmlz']['indx']['grid'] = to_pppack_style(a2_data['nmlz']['indx']['tau'])
    return a2_data


def standarization_ofconc(a2_data):
    """ Saving concentrations in the support PS of support dictionary """
    aux_dic = OrderedDict()
    for i in a2_data:
        evol_tuple = a2_data[i]['conc'].keys()
        if len(evol_tuple) != 1:
            raise RuntimeError('too many tuples for conc')
        evol_tuple = evol_tuple[0]
        aux_dic[i] = a2_data[i]['conc'][evol_tuple]
    return aux_dic


def elimination_ofconc(a2_data):
    """elimination of the concentration from I """
    for data in a2_data.values():
        data.pop('conc')
    return a2_data


def efective_real_state_space(a2_data, num_dim):
    """
    # this can be done in one line, though the changing dic size in interation error has to be dealt with
    # Deleting non samples state parameters
    """
    delete_list = []
    for key in a2_data:
        if len(a2_data[key]) == 1:
            delete_list.append(key)

    for key in delete_list:
        a2_data.pop(key)
    del delete_list

    num_dim_apollo = len(a2_data.keys())
    if num_dim != num_dim_apollo:
        raise RuntimeError('Required dimensions are not present in the A2 file')
    return a2_data


def ordered_real_state_space(a2_data, py_order, a2_order):
    """
    Ordering real-state space according to user defined ordering scheme
    """
    aux_dic = OrderedDict()
    for py_index, key in enumerate(py_order):
        if key in a2_data:
            try:
                a2_index = a2_order.index(key)
            except ValueError:
                a2_index = None
            aux_dic[StateParPickable(key, py_index, a2_index)] = a2_data[key]

    return aux_dic


def normalize_state_space(a2_data):
    """
    Normalization of PS from real values to the hyper-cube [0,1]^d
    """
    aux_dic = OrderedDict()
    for key, vec in a2_data.iteritems():
        _, aux_dic[key] = normalize(vec)
    return aux_dic


def nonphysicalxs_remotion(a2_data,res_nufi_removal):
    """
    Elimination of non-wanted data

    todo=unitest checks
    i=macr doesnt have nufi
    any isotope has higher order anisotropy
    if no inosotropy or no macr no error is produced
    all remaining xs values remain the same
    """
    for i in a2_data['I'].keys():
        if i=='MACR' and res_nufi_removal==True:
            if 'nufi' in a2_data['I'][i]['R'].keys():
                a2_data['I'][i]['R'].pop('nufi')
        for r in a2_data['I'][i]['R'].keys():
            if any(x in r for x in ['111', '112', '122', '212', '222', '211', '322',
                                    '321', '312', '311', '221', '121']):
                a2_data['I'][i]['R'].pop(r)
    return a2_data


def isotopes_remotion(a2_data, retained_isotopes, residual_rebuild=False):
    """
    Elimination of isotopes not in list

    todo=unitest checks
    data only has wanted isotopes
    all remaining xs values remain the same
    """
    if len(retained_isotopes) != len(a2_data['I'].keys()) and residual_rebuild==True:
        a2_data = f_residual_rebuild(a2_data,retained_isotopes)
    for i in a2_data['I'].keys():
        if i not in retained_isotopes:
            a2_data['I'].pop(i)
    if a2_data['I']:
        return a2_data
    else:
        raise RuntimeError('No isotope was retained')


def f_residual_rebuild(a2_data,retained_isotopes,do_rebuild=True):
    """An attempt in rebuilding the proper residual is done. Though sucesfull, the MACR nufi does not really belong to the problem.
        With high specialization this may be a negative value, its excluded from the calculation. Its not trivial task to define when to start considering as the number os specialized isotopes gets reduced
        In addition is also can pose problems in error analysis, specially for relative errors. All this issues wold need to be properly handled before an effective inclsuion of nufi RES
    """
    if do_rebuild!=True:
        return a2_data
    if 'MACR' not in a2_data['I'].keys():
        raise RuntimeError('cannot update MACR as not present')
    if 'MACR' not in retained_isotopes:
        raise RuntimeError('cannot update MACR as not retained, are you sure about this?')
    for i in a2_data['I']:
        if i!='MACR' and i!='MACRT' and i not in retained_isotopes:
            print i
            for r in a2_data['I'][i]['R']:
                if r in a2_data['I']['MACR']['R'].keys(): #avoiding higher anisotropy xs
                    for g in a2_data['I'][i]['R'][r]:
                        for tup in a2_data['I'][i]['R'][r][g]:
                            if len(set(a2_data['I'][i]['conc'].keys()))!=1:
                                raise ValueError('Problems whith the evolution tuple')
                            sigma_perC=np.multiply(a2_data['I'][i]['conc'].values()[0],a2_data['I'][i]['R'][r][g][tup])
                            a2_data['I']['MACR']['R'][r][g][tup]=np.add(a2_data['I']['MACR']['R'][r][g][tup],sigma_perC)
    return a2_data


def reaction_remotion(a2_data, retained_reactions):
    """
    Elimination of reactions not in list
    todo=unitest checks
    data only has wanted isotopes
    all remaining xs values remain the same
    """
    for i in a2_data['I'].keys():
        for r in a2_data['I'][i]['R'].keys():
            if r not in retained_reactions:
                a2_data['I'][i]['R'].pop(r)
    return a2_data


def group_remotion(a2_data, retained):
    """
    Elimination of groups not in list

    todo=unitest checks
    data only has wanted isotopes
    all remaining xs values remain the same
    """
    for i in a2_data['I'].keys():
        for r in a2_data['I'][i]['R'].keys():
            for g in a2_data['I'][i]['R'][r].keys():
                if g not in retained:
                    a2_data['I'][i]['R'][r].pop(g)
    return a2_data


def xs_exists(i, r, g):
    """
    This functions determins tha existance of a given combination of irg. Its used for determining what is consider in the calculation of MXS so any change/error here can have grave consequences difficult to track down

    this function permits to analysise only proper combination of isotope and reaction
    """
    # all istopoes
    act_i = ['U234', 'U235', 'U236', 'U238', 'PU238', 'PU239',
             'PU240', 'PU241', 'PU242', 'NP237', 'AM241', 'AM243']
    fp_i = ['RH103', 'CS133', 'ND143', 'ND145', 'GD155', 'MO95', 'TC99', 'RU101', 'AG107', 'AG109', 'SM147', 'SM149', 'SM150',
            'SM151', 'SM152', 'EU153', 'XE135', 'I135', 'IN115', 'CD106', 'CD108', 'CD110', 'CD111', 'CD112', 'CD113', 'CD114', 'CD116', 'B10', 'B11']

    if i != None and i not in act_i and i not in fp_i and i not in ['MACR', 'MACRT']:
        raise ValueError('Update of iso lists required not present iso', i)

    # Exclusive to cathegory reac
    act_r = ['fiss', 'nufi', 'spec']
    macr_r = ['ener', 'difc', 'tota']

    # for i=None
    if r != None and g != None:
        if 'tran' in r:
            if r[5] == '1' and g != '1':
                return False
            if r[5] == '2' and g != '2':
                return False

    # for i!=None
    if i != None:
        if i in act_i:
            if r in macr_r:
                return False

        if i in fp_i:
            if r in macr_r:
                return False
            if r in act_r:
                return False

        if 'MACR' in i:
            if 'tran2' in r:
                return False
            if 'tran3' in r:
                return False

        # excs are reaction n2n, n3n,... If the iso has high abso, then it doesnt shouw this r
        if i == 'GD155' or i == 'SM150' or i == 'XE135' or i == 'I135' or i == 'XE135' or i == 'B10':
            if r == 'excs':
                return False

    return True


def all_xs(tran_flag=None):
    """
    Standard igr is given here, present in the base irradiation
    """
    xs_ofinterest = {}
    if tran_flag == 'tran':
        xs_ofinterest["r"] = ['abso', 'fiss', 'nufi',
                              'spec', 'tran', 'ener', 'difc', 'tota', 'excs']
    if tran_flag == 'tranXYZ':
        xs_ofinterest["r"] = ['abso', 'fiss', 'nufi', 'spec', 'ener', 'difc', 'tota', 'excs',
                              'tran121', 'tran221', 'tran011', 'tran311', 'tran012', 'tran112', 'tran122', 'tran322', 'tran321', 'tran212', 'tran211', 'tran312', 'tran111', 'tran222', 'tran021', 'tran022']
    xs_ofinterest["g"] = ['1', '2']

    return xs_ofinterest


def do_H5F_2_PY(AX_dic, tag, d):
    """
    routine
    Function that utilizes the class for xs reading. This could eventuallt be proper integrated with the class, thus defining a proper hierarchy of class/subclass
    """

    # accesing xs file
    ps1 = xs_data(AX_dic['path']['file_path'], AX_dic['A2'][tag]['info']['xs_folder'], AX_dic['A2'][tag]['info']['xs_file'],
                  AX_dic['path']['sbr_path'], AX_dic['path']['sbr_file'])  # path for xs and sbr is defines
    ps1.get_phase_space(grid_flag='FG')
    # the auxiliary files are generated with sbr. if generate_out_flag='yes'
    # the *.out files are   generated.
    grid_flag = 'FG'
    ps1.xs_auxiliar_file_generator(AX_dic['A2'][tag]['info']['generate_out_flag'], AX_dic['A2'][tag]['info']['flag_FG2semiFG'], grid_flag,
                                   AX_dic['path']['out_folder'], AX_dic['A2'][tag]['info']['out_alias'])  # grid_flag is required, options; 'SG', 'FG'
    domain_ofinterest = cp.deepcopy(ps1.phase_space)
    xs_ofinterest, domain_ofinterest = domain_reduction(domain_ofinterest, d, AX_dic['A2'][
        tag]['info']['evol_vec'], ps1.order)
    IRG = []
    for key in xs_ofinterest.keys():
        IRG.append('_' + str(len(xs_ofinterest[key])))
    AX_dic['A2'][tag]['info']['IRG'] = ''.join(IRG)
    xs_out, order = ps1.xs_retrival_FG(xs_ofinterest, domain_ofinterest,
                                       AX_dic['path']['out_folder'], AX_dic['A2'][tag]['info']['out_alias'], AX_dic['A2'][tag]['info']['flag_FG2semiFG'])
    conc_dic, fi_dic, k_dic = ps1.cellwise_retrival(domain_ofinterest, AX_dic['path']['out_folder'],
                                                    AX_dic['A2'][tag]['info']['out_alias'], AX_dic['A2'][tag]['info']['flag_FG2semiFG'], AX_dic['A2'][tag]['info']['evol_vec'])

    # The structure of the xs data is here generated
    AX_dic['A2'][tag]['data'] = {}
    AX_dic['A2'][tag]['data']['I'] = xs_out
    AX_dic['A2'][tag]['data']['order_tuple'] = order
    AX_dic['A2'][tag]['data']['PS'] = ps1.domain_ofinterest

    for i in AX_dic['A2'][tag]['data']['I'].keys():
        AX_dic['A2'][tag]['data']['I'][i]['conc'] = conc_dic[i]
    AX_dic['A2'][tag]['data']['fi'] = fi_dic
    AX_dic['A2'][tag]['data']['k'] = k_dic

def do_abso_reconstruction(AX_dic):
    """
    Abso is modified to be abso=abso-excs and excs is removed for good
    """
    for key in AX_dic['A2'].keys():
        for i in AX_dic['A2'][key]['data']['I'].keys():
            if xs_exists(i, 'excs', '1'):
                for g in AX_dic['A2'][key]['data']['I'][i]['R']['abso'].keys():
                    for iota in AX_dic['A2'][key]['data']['I'][i]['R']['abso'][g].keys():
                       # print i,g,iota
                        AX_dic['A2'][key]['data']['I'][i]['R']['abso'][g][iota] = np.array([sg_a - sg_e for sg_a, sg_e in zip(
                            AX_dic['A2'][key]['data']['I'][i]['R']['abso'][g][iota], AX_dic['A2'][key]['data']['I'][i]['R']['excs'][g][iota])])
                AX_dic['A2'][key]['data']['I'][i]['R'].pop('excs')


def domain_reduction(domain_ofinterest, d, evol_tuple, evol_order_tuple):
    if evol_order_tuple[1] != 'DENSITE_MODERATEUR':
        raise ValueError(
            'Im not projecting in an evol value for the MD, the value Im using may not even exist')
    # this should be the standard
    #"""
    if d < 4:
        domain_ofinterest['DENSITE_MODERATEUR'] = [domain_ofinterest['DENSITE_MODERATEUR'][
            evol_tuple[0][1] - 1]]  # A2 indexation starts at 1 and PY at 0
    if d < 3:
        domain_ofinterest['CONCENTRATION_EN_BORE'] = [
            domain_ofinterest['CONCENTRATION_EN_BORE'][evol_tuple[0][2] - 1]]
    if d < 2:
        domain_ofinterest['TEMPERATURE_COMBUSTIBLE'] = [
            domain_ofinterest['TEMPERATURE_COMBUSTIBLE'][evol_tuple[0][0] - 1]]
    #"""
    xs_ofinterest = {}

    xs_ofinterest["i"] = ['U234', 'U235', 'U236', 'U238', 'PU238', 'PU239', 'PU240',
                          'PU241', 'PU242', 'NP237', 'AM241', 'AM243', 'RH103', 'CS133', 'ND143', 'ND145',
                          'GD155', 'MO95', 'TC99', 'RU101', 'AG109', 'SM147', 'SM149', 'SM150',
                          'SM151', 'SM152', 'EU153', 'XE135', 'I135', 'MACR', 'MACRT']  # ,'MACRT'
    xs_ofinterest["r"] = ['abso', 'nufi', 'tran', 'excs']
    xs_ofinterest["g"] = ['1', '2']
    return xs_ofinterest, domain_ofinterest


class xs_data(object):
    r""" xs class for reading few-group xs (A2 outpu). The class functionalities are:

    * Recive input paths (xs and sbr file) and output path (for generation of auxiliary .out file and logs)
    * Obtain the phase-space from xs file (Names, discretization and values)
    * Retrive xs data from xs file. This is done by using the sbr.exe and generating auxiliary files (*.out). After auxiliary files generation xs of interest are retrived and stored in an orderd way by matching the points of calculation with the corresponding phase space values
    * 3 phase space shcemes are considered. Full grid (everything is proper cartesian product), semi-FG (some part pof the cartesian product are replaced with user imposed functions) and SG (boils down to node-wise discretization: 1D)

    Note: At the moment of writing it is necesary to check the validity of all section where 'USER IMPOSED' tag can be read.
    """

    def __init__(self, file_path, xs_folder, xs_file, sbr_path, sbr_file):
        "Names and paths for the xs (h5f) and sbr file is recived"

        xs_path = os.path.join(file_path, xs_folder, xs_file)
        sbr_path = os.path.join(sbr_path, sbr_file)

        if not os.path.isfile(xs_path):
            raise ValueError("could not find the xs file in %s" % xs_path)
        if not os.path.isfile(sbr_path):
            raise ValueError("could not find the sbr file in %s" % sbr_path)

        self.working_path = os.getcwd() + '/'
        self.file_path = file_path  # path to input file
        self.xs_folder = '/' + xs_folder + '/'  # path to input file
        self.xs_file = xs_file  # path to input file
        self.xs_path = xs_path  # includes the name of xs file, i.e. /home/user/.../xs_file
        self.sbr_file = sbr_file  # Name of file given by user
        self.sbr_path = sbr_path  # includes the name of sbr file, i.e. /home/user/.../sbr_file

    def get_phase_space(self, grid_flag):
        r"""
        Xs's phase-space is obtained directly from the h5f and store in the class. This permits a general tratment for the full range of the phase-space without user intervention.

        Self.
        * N      : a vector containing discretization (cardinality) for all the dimensions di
        * d      : number of considered dimensions
        * NPAR   : number of parameters in the lattice files
        * phase_space: dictionary where the phase space is stored preserving the original
                       keys from the lattice file for 1\leq d_i\leq d
        * order: saves the order of the different keys in the lattice file
        """

        f = h5py.File(self.xs_path, 'r')
        self.N = f['paramdescrip']['NVALUE'].value  # det maximum range Ni for each d_i
        phase_space = {}
        order = {}
        NPAR = f['paramdescrip']['NPAR'].value[0]
        for di in range(NPAR - 1):
            di_name = f['paramdescrip']['PARNAM'].value[di]  # get names for dimensions. Starts at 0
            # get values for dimensions. Starts at 1. e.g. 'BURNUP': array([  0., 9.35253143, 18.70503998,..
            # Is saved as a np.array, of floats64 FORTRAN-contiguous
            phase_space[di_name] = np.array([float(val) for val in f['paramvaleurs'][
                'pval       %d' % (di + 1)].value], order='F')
            order[di] = di_name  # e.g. '7': 'BURNUP'

        iso_aux = []
        # just concatenate those two
        for iso in f['contenu']['NOMISO'].value[:]:
            iso_aux.append(iso)
        for iso in f['contenu']['NOMMAC'].value[:]:
            iso_aux.append(iso)
        f.close()
        self.iso_A2 = iso_aux

        # USER IMPOSED: Non-independant variables set to [0].
        """
        *Do not eliminate them, this will bring problems with the cartesin product later one
        *if instead of '[phase_space['PHASE'][0]]' (which is equal to 1) just '[1]' is written then np.where() does not recognize the value.

        This two problems rise from the decision of defining the 'space of interest' as a subset from the 'phase space' which in time is read directly from the H5F file. Later several comparisons are made between the two. The upside is the need for no explicit declaration of the phase-space thus minimizing chances of un-noticed error in domain assignation.
        """
        if 'PHASE' in phase_space.keys():
            phase_space['PHASE'] = [phase_space['PHASE'][0]]
        if 'BURNUPstep' in phase_space.keys():
            phase_space['BURNUPstep'] = [phase_space['BURNUPstep'][0]]

        if grid_flag == 'SG':  # major update required
            """
            In contras to FG, the stored values in the concatenated SAPHYB file only considers different burnup steps, i.e a set of values [0, 500, 500, 100] are stored as [0, 500, 100]. Two posibilities remain, read the BURNUP value from the single XS files separatly or load a pickeled object with the phase space. The second option was implemented.
            """
            with open(self.file_path + self.xs_folder + 'phase_space.pickle', 'rb') as handle:
                phase_space_pk = pickle.load(handle)
                phase_space_pk.pop('a')
                phase_space_pk.pop('d')
                phase_space_pk.pop('l')
                phase_space_pk.pop('BURNUP_evol')
                phase_space_pk.pop('BURNUP_steps')
            phase_space = phase_space_pk

        self.phase_space, self.order, self.d, self.NPAR = phase_space, order, len(order), NPAR

    def __str__(self):
        print "structure of tree and content"

    def xs_auxiliar_file_generator(self, generate_out_flag, flag_FG2semiFG, grid_flag, out_folder, out_alias):
        r"""
        For extracting xs a phase space grid must be specified. For a full grid or semi-full grid(FG) composed by a cartesian product of the discretized domain d_i 1\leq d_i \leq d then FG='yes'.
        If another type of grid has been considered (e.g. sparse grid) then an input table of valid points need to be additionaly provided to generate sbr.dat, run sbr and obtain the auxiliary *.out and FG='no'
        generate_out_flag: determins if the calculation of the auxiliary files will be carried ou
        Note: at moment of writing passing trough 'xs_auxiliar_file_generator' is mandatory even if auxiliary files are not generated. Since in this step the conversion table is built where the relation between the nodes indexes and the .out files gets fixed.

        Note: If auxiliary files already exist and generate_out_flag='yes' then they are deleted. Generating .out files can be time consuming.
        """
        if grid_flag == 'FG':
            self._FG_saphyr(generate_out_flag, flag_FG2semiFG, out_folder, out_alias)

        elif grid_flag == 'SG':
            self._SG_saphyr(generate_out_flag)
        else:
            print "Method for getting xs from non FG not imlemented yet. Input file of valid state point is needed"

    def _FG_saphyr(self, generate_out_flag, flag_FG2semiFG, out_folder, out_alias):
        r"""xs auxiliary files generation for full grid:=points of calculations are the cartesian product of the discretized domain. Also semi-full grid can be consider with USER IMPOSED rules

        conversion table: permits to link the index-point of calculation with the index-touple of the disretizied phase-space. list of pairs ((phase-space index set),#.out file), e.g. ((1, 1, 1, 1, 1, 1, 1, 15), 15). Permits to link the auxiliary files with the phase-space. Due to this step running 'xs_auxiliar_file_generator' is mandatory

        Note: At the moment of writing the index of the sequence of the tuple in the conversion table (i.e. ..., (1, 1, 1, 1, 1, 1, 1, 15), (1, 1, 1, 1, 1, 1, 1, 16),...) is enough for indexating .out files. However this has been left since later a custum indexation could prove usefull (for very large files maybe).
        """
        # Generate folder for *.out files
        os.chdir(out_folder)
        # consistacy of folder declaration with content of sbr.dat file is needed
        if not os.path.isdir(out_folder + out_alias):
            print "making dir " + out_folder + out_alias
            os.mkdir(out_folder + out_alias)
        else:
            if generate_out_flag == 'YES':
                # print "rm "+out_folder+out_alias+"/*.out"+"> /dev/null 2>&1"
                os.system("rm " + out_folder + out_alias + "/*.out" + "> /dev/null 2>&1")

        # A general FG is generated
        FG = []
        for di in range(self.d):
            FG.append(range(1, self.N[di] + 1))

        # The general FG is adapted to the actual semiFG present in the saphyb files, define by the user
        # The degree of freedom will be reduced to those utilize by A2 user in the reprise
        FG = self.FG2semiFG(FG, flag_FG2semiFG)
        # print FG; sys.exit()

        # Generate valid touples as cartesian product of lists contain in FG:
        # tuple=\times_1^{d}\{x_j\text{ }1\leq j\leq N_i\}
        conversion_table = []
        point_calc = 1

        # The sbr.exe needs to be in the same folder as the input file. IDK why...
        if generate_out_flag == 'YES':
            os.system("cp " + self.sbr_path + " " + out_folder + out_alias)
            print 'sbr copied'
        if generate_out_flag == 'YES':
            print "time for generating all out files about [horus, days]", (40. / (60. * 60.)) * len([p for p in itertools.product(*FG)]), (40. / (60. * 60. * 24)) * len([p for p in itertools.product(*FG)])

        for tuple_calc in itertools.product(*FG):

            tuple_calc = np.array(tuple_calc)  # required for tuple manipulation
            # With the proper semiFG, the adequate tuple for retriving data is
            # re-build by aconsidering the actual variables used by the user in A2
            tuple_calc = self.tupleFG2tuple_semiFG(tuple_calc, flag_FG2semiFG)
            # Here the relation of the point of calculation in the .out file and the actual tuple in the spahyb is saved
            # Modifications to the tuple_calc needs to be done before this point or
            # un-existant calculation point will be consider later on, producing an
            # error
            conversion_table.append((tuple_calc, point_calc))
            # sbr.exe sbr.dat manipulation, for generating .outfiles
            if generate_out_flag == 'YES':
                print "processing calculation point #", point_calc, " for tuple ", tuple_calc
                # Generation of input for saphyb.exe
                sbr_input = open(out_folder + out_alias + '/' + 'sbr.dat', 'w')
                sbr_input.write(self.xs_path + "\n")
                sbr_input.write(out_folder + out_alias + "/out.lst" + "\n")
                sbr_input.write(out_folder + out_alias + "/" + str(point_calc) + ".out" + "\n")
                sbr_input.write(str(self.NPAR) + "\n")
                # Formating tuple_calc. e.g from '(1, 1, 1, 1, 1, 1, 1, 13)' to '1 1 1 1 1 1 1 13'
                formated_tuple_calc = ""
                for pi in tuple_calc:
                    formated_tuple_calc = formated_tuple_calc + str(pi) + " "
                sbr_input.write(formated_tuple_calc + "\n")
                sbr_input.close()
                os.chdir(out_folder + out_alias + "/")
                p = subprocess.Popen('./sbr.exe sbr.dat>>sbr.out', shell=True,
                                     stdin=subprocess.PIPE, universal_newlines=True)
                # 0.3 for file of 2600 points max (at 3900 crashes)
                # 1.0 achives equilibrium for a file of 3900 points. (i.e. multiple process running without overlaping), less could be considered
                # 4.0 achives equilibrium for a file of 8000 points. (i.e. multiple process running without overlaping), less could be considered
                # I try 20.0 for a file of 48000 points
                # time of execution of one calc point for one saphyb of 3900 is 17 sec
                # time of execution of one calc point for one saphyb of 7884 is 22 sec
                # time of exec of one calc point of 50000 is 109 sec
                os.system("sleep 20.0")
            point_calc += 1

        # Cheking that non-existing points have not been requested (if sbr.out file exists)
        if generate_out_flag == 'YES':
            os.system("sleep 500.0")  # give time for the last points to be obtained
        if os.path.isfile(out_folder + out_alias + '/' + "sbr.out"):
            fout = open(out_folder + out_alias + '/' + "sbr.out", 'r')
            for lines in fout:
                if lines.find('non existing state point having parameters') != -1:
                    raise ValueError(
                        'requested non-existent phase-space point. I suggest taker a deeper close to functions FG2semiFG and  tupleFG2tuple_semiFG %s' % lines)

        os.chdir("..")

        if generate_out_flag == 'YES':
            fout = open(out_folder + out_alias + '/' + "tuple_calc.out", 'w')
            fout.write('[Calcul tuple,.out file')
            for tuple_calc in conversion_table:
                fout.write(str(tuple_calc))
                fout.write('\n')
            fout.close()

        self.conversion_table = conversion_table

    def FG2semiFG(self, FG, flag_FG2semiFG):
        """
        Here the user imposed rules for the saphyb files are explicitly commented/uncommented. In the future a more elaborated system can be considered.

        As a general rule the idea is to reduced the proposed FG generated by the  cartesian product of what A2 consider the whole phase space domain to the actual minimum degrees of freedom, possibly mimiquing the A2 strategy for generating the branch calculation.
        Later other function will reconstrcut the proper tuple considering the parameters allow to vary here.
        """
        # USER IMPOSED
        # Actual FG
        if 'JUSTE FG' in flag_FG2semiFG:
            return FG
        # Antonio's files where a 'phase' parameter is used to permit changin
        # evolution condition (CR insertion) at 2 points: 15,30 MWd/tU thus
        # defining 3cicles
        if 'ESTEBAN' in flag_FG2semiFG:
            """
            Esteban's files. 'phase' and 'burnup step' and 'burnup' depends on a single parameter. No additional point required
            Example
            input FG:[[1], [1], [1], [1], [1], [1], [1], [1], [1, 2,..., 70]]
            Output FG:[[1], [1], [1], [1], [1], [1], [1], [1], array([ 1,  2,..., 70, 71, 72])]
            """
            FG[6] = [1]  # the PHASE is a dependant variables and thus does not participates in the cartesian product
            FG[7] = [1]  # the BURNUP_step is a dependant variable...
        # Antonio's files where a 'phase' parameter is used to permit changin
        # evolution condition (CR insertion) at 2 points: 15,30 MWd/tU thus
        # defining 3cicles
        if 'ANTONIO' in flag_FG2semiFG:
            """
            Antonio's files. 'phase' and 'burnup step' and 'burnup' depends on rules over a special index
            Example
            input FG:[[1], [1], [1], [1], [1], [1], [1], [1], [1, 2,..., 70]]
            Output FG:[[1], [1], [1], [1], [1], [1], [1], [1], array([ 1,  2,..., 70, 71, 72])]
            Two points are added, as A2 calculation at 0 burnup are mandatory, and at each new cicle a new 0 burnup calcvulation takes place
            """
            FG[6] = [1]  # the PHASE is a dependant variables and thus does not participates in the cartesian product
            FG[7] = [1]  # the BURNUP_step is a dependant variable...
            # the phase depended on the index in burnup by this relationship. The number 24,48 are the indexed defined by the discretization criterion of burnup chosen by the user of A2
            # print FG
            if FG[8][-1] <= 24:
                ph = 1
            if FG[8][-1] > 24 and FG[8][-1] <= 48:
                ph = 2
            if FG[8][-1] > 48:
                ph = 3
            # print ph
            # it adds indexes as required by ph. if ph=1 then nothing gets added, if
            # ph=3 and FG[8][-1]=70 then '71' and '72' get added
            FG[8] = np.append(FG[8], range(FG[8][-1] + 1, FG[8][-1] + ph))
            if 'CR' in flag_FG2semiFG:
                """
                In reference calculation where CR position may change, the lecture of the phase space given a CR value of [1,2] i.e. in and out. However no branch calculation acctualy takes place, and un-existing points are requested.
                """
                FG[4] = [1]
            if 'NO BURNUP' in flag_FG2semiFG:
                FG[8] = [1]
        # Format: Making sure that index are int
        return [[int(value) for value in vec] for vec in FG]

    def tupleFG2tuple_semiFG(self, tuple_calc, flag_FG2semiFG):
        """
        The tuple reconstruction for reading saphyb files from the actual degrees of freedom permitted by FG2semiFG
        """
        # USER IMPOSED
        """
        Full grid
        """
        if 'just FG' in flag_FG2semiFG:
            """
            Example:
            Input:  [ 1  1  1  1  1  1  1  1 72]
            Output: [ 1  1  1  1  1  1  1 70 70 70]
            """
            tuple_calc = np.append(tuple_calc, [tuple_calc[-1], tuple_calc[-1]])

        if 'ESTEBAN' in flag_FG2semiFG:
            """
            Example:
            Input:  [ 1  1  1  1  1  1  1  1 72]
            Output: [ 1  1  1  1  1  1  1 1 70 70]
            last 3 indexes are phase=1 burnup=burnupstep
            """

            tuple_calc[7] = tuple_calc[-1]
            tuple_calc = np.append(tuple_calc, [tuple_calc[-1]])
        if 'ANTONIO' in flag_FG2semiFG:
            """
            Example:
            Input:  [ 1  1  1  1  1  1  1  1 72]
            Output: [ 1  1  1  1  1  1  3 24 70 70]
            """
            aux = tuple_calc[8]
            if tuple_calc[8] <= 24:
                ph = 1
            if tuple_calc[8] > 24 and tuple_calc[8] <= 48:
                ph = 2
            if tuple_calc[8] > 48:
                ph = 3
            tuple_calc[6] = ph
            tuple_calc[7] = aux - 24 * (ph - 1)
            tuple_calc[8] = aux - (ph - 1)
            tuple_calc = np.append(tuple_calc, [aux - (ph - 1)])
            # this is for the reference calculation of CR with no branch with partial
            # insertion during the cicles
            if 'CR' in flag_FG2semiFG and 'FIRST CYCLE' in flag_FG2semiFG:
                if ph == 1:
                    tuple_calc[4] = 2

            if 'CR' in flag_FG2semiFG and 'LASR CYCLE' in flag_FG2semiFG:
                if ph == 3:
                    tuple_calc[4] = 2

            if 'CR' in flag_FG2semiFG and 'MIDDLE CYCLE' in flag_FG2semiFG:
                if ph == 2:
                    tuple_calc[4] = 2

        return tuple_calc

    def xs_retrival_FG(self, xs_ofinterest, domain_ofinterest, out_folder, out_alias, flag_FG2semiFG):
        r""" For the xs of interest (isotope,reac,group) and the domain of interest in the phase space, the required calculation points are read and xs builded as objects in a xs dictionary. Additionally if normalization_flag='yes' the domain is scales to [0,1] for all dimensions. Only the files required for the domain_ofinterest are read.
        Expected input:
        xs_ofinterest  :  {'i': ['U235'], 'r': ['abso'], 'g': ['2']}
        """
        self.iso_read = xs_ofinterest['i']
        # only isotopes of interest are going to be read. However, iso_A3 and
        # iso_read should be the same if macroscopic XS are going to be
        # calculated.
        # A list is generated. Each element is another list with the index
        # positions of the requested domain in the phase space. e.g. [[3], [1],
        # [1], [1], [1], [1], [1], [1, 2, 3, 4, 5]]. Self.order establishes the
        # link between the phase space index of a given dimension and its names
        # (keys). Any manipulation on the domain of interest must not invalidate
        # the search np.where(), otherwise empty arrays (array()) may come up.
        idx_tuple_calc = []
        for di in range(self.d):
            idx_tuple_calc.append([np.where(val == self.phase_space[self.order[di]])[0][
                                  0] + 1 for val in domain_ofinterest[self.order[di]]])
        # print idx_tuple_calc
        idx_tuple_calc = self.FG2semiFG(idx_tuple_calc, flag_FG2semiFG)
        # print idx_tuple_calc;sys.exit()
        order = [self.order[i] for i in range(0, 6)]
        # I want to locate XS for index in phase space. So a USER DEFINED set of indexes is considerd
        # The parametrization is on iota, so only [0:6] is considered, but I do need to apply the rules on FG and tupleFG2tuple_semiFG for assuring consistancy of variables.
        #'''
        # generating anisotropy vector
        # This can be passed further up if in the future many files have different
        # number of groups or anisotropy levels
        anysotropy = 3
        anysotropy_vec = [str(lvl) for lvl in range(anysotropy + 1)]
        groups = 2
        groups_vec = [str(lvl) for lvl in range(1, groups + 1)]
        # generation of xs dictionary
        xs_dic = {}
        for i in xs_ofinterest['i']:
            xs_dic[i] = {}
            xs_dic[i]['R'] = {}
            for r in xs_ofinterest['r']:
                if r != 'tran':
                    if xs_exists(i, r, None):
                        xs_dic[i]['R'][r] = {}
                        for g in xs_ofinterest['g']:
                            xs_dic[i]['R'][r][g] = {}
                            for tuple_i in itertools.product(*idx_tuple_calc):
                                aux = tuple(self.tupleFG2tuple_semiFG(
                                    np.array(tuple_i), flag_FG2semiFG))
                                # print aux
                                xs_dic[i]['R'][r][g][aux[0:6]] = []
                else:
                    """
                    tran XS are saved indexed as 'tran'+'anisotropy level'+'input group'+'output group'
                    level 0 are the standard scaterring xs for a whole assembly flux. So:
                    tran011=\sigma_{1->1},tran012=\sigma_{1->2},tran021=\sigma_{2->1},tran022=\sigma_{2->2}

                    Note: scaterring xs for iso=MACR and anisotropy>1 is generated, i.e. tran2** and tran3** but then they won't be filled with anything
                    """
                    for p in anysotropy_vec:
                        for g1 in groups_vec:
                            for g2 in groups_vec:
                                # print r+p+g1+g2
                                if xs_exists(i, r + p + g1 + g2, None):
                                    xs_dic[i]['R'][r + p + g1 + g2] = {}
                                    xs_dic[i]['R'][r + p + g1 + g2][g1] = {}
                                    for tuple_i in itertools.product(*idx_tuple_calc):
                                        aux = tuple(self.tupleFG2tuple_semiFG(
                                            np.array(tuple_i), flag_FG2semiFG))
                                        xs_dic[i]['R'][r + p + g1 + g2][g1][aux[0:6]] = []
        # From the list of required indices of d dimensions a list of tuples is
        # build. For the requested tuples, a point of calculation in the auxiliary
        # *.out files is found by self.conversion_table. The condition for a
        # requesting a point of calculation is a match between the tuple and the
        # available touples in the phase space. e.g [49, 50, 51, 52, 53, 54, 55,
        # 56, 57, 58]. If user-imposed specific 'non FG' whese consider in the
        # conversion table generation here they need to be considered as well

        point_calc = None
        for tuple_i in itertools.product(*idx_tuple_calc):
            # USER IMPOSED: the conversion table saves user defined relation in the
            # indexes of the nodes
            tuple_i = self.tupleFG2tuple_semiFG(np.array(tuple_i), flag_FG2semiFG)
            # print tuple_i
            # for the requested tuple_i the corresponding .out file is found
            for i in range(len(self.conversion_table)):
                if all(tuple_i == self.conversion_table[i][0]):
                    # the conversion table permits to consider custom naming of .out files
                    point_calc = self.conversion_table[i][1]
                    break  # calculation points are unique. After the first match the search for that tuple is abandoned
            if i == len(self.conversion_table):
                raise ValueError(
                    'a point not existing in the .out files has been requested. tuple=', tuple_i)

            # un-comment for locating specific .out files in the xs reading process
            """
            if all(tuple_i==[ 2,  2,  1,  1,  2,  1 , 1, 1, 1, 1]) or all(tuple_i==[ 2,  2,  1,  1 , 2 , 1 , 1 ,24 ,24, 24]):
                print tuple_i, point_calc

            if all(tuple_i==[ 2,  2,  1,  1 , 1 , 1 , 2 ,1 ,24, 24]):
                print tuple_i, point_calc
            """

            # Access auxiliary *.out files
            fout = open(out_folder + out_alias + '/' + str(point_calc) + ".out", 'r')
            iso = None

            for line in fout:
                # Detect isotopes specification
                if line.find('isotope') != -1:
                    iso = line.split()[1]
                    tran_counter = 0

                # Append only xs of interest. tran is a special case and treated as group independent
                # print xs_ofinterest;sys.exit()
                if iso in xs_ofinterest["i"]:
                    for reac in ['abso', 'fiss', 'nufi', 'spec', 'tran', 'ener', 'difc', 'tota', 'excs']:
                        # A xs may not be present, this automaticly handled by line.find(reac)!=-1
                        # A xs may be present but not wanted, this is handled by: reac in xs_ofinterest["r"]
                        # A xs may be unphysical (nufi in MACR) this is handle by
                        # xs_exists(iso,r,None)
                        if line.find(reac) != -1 and reac in xs_ofinterest["r"] and xs_exists(iso, reac, None):
                            if reac != 'tran':
                                # print iso, reac,xs_dic[iso]['R'].keys(),  xs_exists(iso,reac,None)
                                if '1' in str(xs_ofinterest["g"]):
                                    xs_dic[iso]['R'][reac]['1'][
                                        tuple(tuple_i[0:6])].append(float(line.split()[1]))
                                if '2' in str(xs_ofinterest["g"]):
                                    xs_dic[iso]['R'][reac]['2'][
                                        tuple(tuple_i[0:6])].append(float(line.split()[2]))
                            else:
                                # this is for P3 anisotropy. Associating a group preservs structure
                                # of dictionary.
                                xs_dic[iso]['R'][
                                    reac + str(tran_counter) + '1' + '1']['1'][tuple(tuple_i[0:6])].append(float(line.split()[1]))
                                xs_dic[iso]['R'][
                                    reac + str(tran_counter) + '1' + '2']['1'][tuple(tuple_i[0:6])].append(float(line.split()[3]))
                                xs_dic[iso]['R'][
                                    reac + str(tran_counter) + '2' + '1']['2'][tuple(tuple_i[0:6])].append(float(line.split()[2]))
                                xs_dic[iso]['R'][
                                    reac + str(tran_counter) + '2' + '2']['2'][tuple(tuple_i[0:6])].append(float(line.split()[4]))
                                tran_counter += 1
            fout.close()
        self.domain_ofinterest = domain_ofinterest
        for i in xs_dic.keys():
            for r in xs_dic[i]['R'].keys():
                for g in xs_dic[i]['R'][r].keys():
                    for iota in xs_dic[i]['R'][r][g].keys():
                        if len(xs_dic[i]['R'][r][g][iota]) != len(domain_ofinterest['BURNUP']):
                            print i, r, g, iota
                            raise ValueError("empty entries for")

        # if zero values are prefared to inexistent data (for isotopes associated
        # to CR and things like that)
        AD_HOC_ZERO = 'no'
        i0 = xs_dic.keys()[0]
        r0 = xs_dic[i0]['R'].keys()[0]
        g0 = xs_dic[i0]['R'][r0].keys()[0]
        iota0 = xs_dic[i0]['R'][r0][g0].keys()[0]
        aux = len(xs_dic[i0]['R'][r0][g0][iota0])

        if AD_HOC_ZERO == 'yes':
            for i in xs_dic.keys():
                for r in xs_dic[i]['R'].keys():
                    for g in xs_dic[i]['R'][r].keys():
                        for iota in xs_dic[i]['R'][r][g].keys():
                            print iota, len(xs_dic[i]['R'][r][g][iota])
                            if len(xs_dic[i]['R'][r][g][iota]) == 0:
                                xs_dic[i]['R'][r][g][iota] = np.zeros(aux)

        return xs_dic, order

    def cellwise_retrival(self, domain_ofinterest, out_folder, out_alias, flag_FG2semiFG, evol_vec):
        """
        Permits the extraction of cell-wise, burnup-wise data for the domein of interest:
        * Pu concentration
        * Fast and thermal flux


        Here there is a problem related to the suple on which the concentrationss are searched for. The evol vector is used:
            *First of all if the PS of interest does not include de evol vector the conc dic will be empty
            *Second, it is not acctualy the evol tuple that is intresting but only a few variables such as CR that may be of interest. At the present moment I see without discrimination the first four values fo the tuple, that many times have no sense whatsoever as they are fuel temperature boron concentration and so on... To correct this variables that had an impact on the concentration should be further specifice in the diagram of the information, and is not worth it, at the moment.

        Note that the conocenctration is obtained only for the evol vector which must be included in the tuple being extracted
        """
        idx_tuple_calc = []
        for di in range(self.d):
            idx_tuple_calc.append([np.where(val == self.phase_space[self.order[di]])[0][0] + 1
                                   for val in domain_ofinterest[self.order[di]]])  # +1 since numeration of .out starts at A
        # concentration dictionary
        conc_dic = {}
        conc_dic_aux = {}
        aux = all_xs()
        for key in self.iso_read:
            conc_dic[key] = {}
        # K_inf dictionary
        k_dic = {}
        k_dic['k_inf'] = {}
        k_dic['k_eff'] = {}
        k_dic['b2'] = {}
        # generation of xs dictionary
        fi_dic = {}
        fi_dic['1'] = {}
        fi_dic['2'] = {}
        fi_dic['si'] = {}
        # fluxes for every point in the phase space. So a USER DEFINED set of indexes is considerd
        for tuple_i in itertools.product(*idx_tuple_calc[0:6]):
            fi_dic['1'][tuple_i] = []
            fi_dic['2'][tuple_i] = []
            k_dic['k_inf'][tuple_i] = []
            k_dic['k_eff'][tuple_i] = []
            k_dic['b2'][tuple_i] = []
            fi_dic['si'][tuple_i] = []
            for i in conc_dic.keys():
                if tuple_i[0:4] == evol_vec[0][0:4]:
                    conc_dic[i][tuple_i] = []

        idx_tuple_calc = self.FG2semiFG(idx_tuple_calc, flag_FG2semiFG)
        point_calc = None


        for tuple_i in itertools.product(*idx_tuple_calc):
            # USER IMPOSED: the conversion table saves user defined relation in the
            # indexes of the nodes
            tuple_i = np.array(tuple_i)
            # print tuple_i,self.tupleFG2tuple_semiFG(tuple_i,flag_FG2semiFG)
            tuple_i = self.tupleFG2tuple_semiFG(tuple_i, flag_FG2semiFG)
            # for the requested tuple_i the corresponding .out file is found
            for i in range(len(self.conversion_table)):
                if all(tuple_i == self.conversion_table[i][0]):
                    # the conversion table permits to consider custom naming of .out files
                    point_calc = self.conversion_table[i][1]
                    break  # calculation points are unique. After the first match the search for that tuple is abandoned
            if i == len(self.conversion_table):
                raise ValueError(
                    'a point not existing in the .out files has been requested. tuple=', tuple_i)
            fout = open(out_folder + out_alias + "/" + str(point_calc) + ".out", 'r')

            for line in fout:
                # Saving flux
                if line.find('fluxs') != -1:
                    fi_dic['1'][tuple(tuple_i[0:6])].append(float(line.split()[1]))
                    fi_dic['2'][tuple(tuple_i[0:6])].append(float(line.split()[2]))

                # Saving Kinf
                if line.find('kinf') != -1:
                    k_dic['k_inf'][tuple(tuple_i[0:6])].append(float(line.split()[1]))

                # Saving Keff
                if line.find('keff') != -1:
                    k_dic['k_eff'][tuple(tuple_i[0:6])].append(float(line.split()[1]))

                # Saving b2
                if line.find('b2') != -1:
                    k_dic['b2'][tuple(tuple_i[0:6])].append(float(line.split()[1]))

                # Saving si
                if line.find('si') != -1:
                    fi_dic['si'][tuple(tuple_i[0:6])].append(float(line.split()[1]))

                # Saving conc
                # if the isotope is in the read isotopes
                if line.find('concentration') != -1 and line.split()[1] in self.iso_read:
                    # print tuple_i,evol_vec

                    if str(tuple(tuple_i[0:4])) == str(evol_vec[0][0:4]):
                        conc_dic[line.split()[1]][tuple(tuple_i[0:6])
                                                  ].append(float(line.split()[4]))

            fout.close()

        # formating
        for i in conc_dic.keys():
            for iota in conc_dic[i].keys():
                conc_dic[i][iota] = np.array(conc_dic[i][iota])

        return conc_dic, fi_dic, k_dic


def to_pppack_style(tau):
    # saving the list in the given order
    """
    [(0.0, 0.0, 0.0), 0]
    [(0.00020842045680406521, 0.0, 0.0), 1]
    [(0.00041684053106259416, 0.0, 0.0), 2]
    [(0.00083367698163946816, 0.0, 0.0), 3]
    ...
    Critical  function
    """

    Ni = [len(t) for t in tau]
    aux = []
    # generating grid points with the required order for pppack
    for x, i in zip(itertools.product(*tau), itertools.product(*[range(nx) for nx in Ni])):
        aux.append([x, getidx(i, Ni)])
    grid = []
    for point in sorted([v for v in aux], key=lambda tup: (tup[1])):
        grid.append(point[0])
    return grid
