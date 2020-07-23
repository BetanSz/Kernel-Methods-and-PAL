# --*-- coding:utf-8 --*--

"""
User's input definitions.

"""
from packages.tools_SG_common import *


def user_case():
    """
    case data run
    """
    usedic=OrderedDict()
    usedic['IS_size'] = 'MEDIUM'
    #usedic['IS_size'] = 'BIG'
    usedic['num_dim'] = 2
    usedic['density_regulizer']=0 #extraction rate of test points at to flattern the test histogram. In the thesis 0 and 6
    usedic['A2_case'] = 'Bench_PWR_UO2_HFP'
    usedic['OS_size'] = '_X1I5R4G1_' #smallest
    #usedic['OS_size'] = '_X1I31R6G2_' # full case
    #usedic['OS_size'] = '_X1I30R4G2_' #article case
    usedic['run_case'] = usedic['A2_case'] +usedic['OS_size'] + usedic['IS_size']+'_' + 'd' + str(usedic['num_dim'])
    #Results
    usedic['less_sermapoints']=None

    # the divisions of each part
    usedic['std']={} #standarization
    usedic['app']={} #calculation of the approximation
    usedic['rlt']={} #evaluation and statistics
    usedic['plt']={} # plotting

    raw_f = '/raw_xs_data/'
    std_f = '/std_xs_data/'
    app_f = '/app_xs_data/'
    rlt_f = '/rlt_xs_data/'

    enviroment_flag='es251151'
    #enviroment_flag='Karim'
    # enviroment_flag='HPC'
    enviroment_flag='home'
    if enviroment_flag=='es251151':
        root_I = '/data/tmplca/eszames/xs_analysis/'#'/volatile/SG_results/'
        root_O = '/data/tmplca/eszames/xs_analysis/'
    if enviroment_flag=='home':
        root_I = '/home/betan/Codes/EstebanCode1.0/Esteban_XS_pack/std/'
        root_O = '/home/betan/Codes/EstebanCode1.0/Esteban_XS_pack/std/'
    if enviroment_flag=='Karim':
        root_I = '/data/tmplca/eszames/xs_analysis/'
        root_O = '/data/tmplca/eszames/xs_analysis/'

    usedic['std']['root_I']=root_I
    usedic['std']['root_O']=root_O
    usedic['std']['f_I']=raw_f
    usedic['app']['root_I']=root_I
    usedic['app']['root_O']=root_O
    usedic['app']['f_I']=std_f
    usedic['rlt']['root_I']=root_I
    usedic['rlt']['root_O']=root_O
    usedic['rlt']['f_I']=app_f
    usedic['plt']['root_I']=root_I
    usedic['plt']['root_O']=root_O
    usedic['plt']['f_I']=rlt_f
    usedic['plt']['f_O']='/opt_xs_data/'

    return usedic

def get_user_TT_division(IS_size,num_dim):

    unit='per'
    #unit='len(idx)'
    if unit =='per':
        """value is % of training from available data, the rest is test"""
        if IS_size == 'MEDIUM':
            if num_dim == 1:
                ret_list = [70]
            if num_dim == 2:
                ret_list = [70, 100]
            if num_dim == 3:
                ret_list = [70, 100, 100]

        if IS_size == 'BIG':
            if num_dim == 1:
                ret_list = [70]
            if num_dim == 2:
                ret_list = [50, 50]
                ret_list = [60, 80]
            if num_dim == 3:
                ret_list = [50, 50, 50] # this split is 5500 points
                ret_list = [60, 50, 50] # BRUTAL AL
                ret_list = [30, 45, 45] # this split still uses all the serma point but is 2000 points in total is the one used for most of the work

    return {unit: ret_list}

def get_user_supp_input(IS_size,num_dim):
    """Typical user data for support budget
    the presence of a serma discret is controlled directly by the flag serma_supp
    """
    unit='per'
    #unit='len(idx)'
    if unit=='len(idx)':
        if num_dim == 1:
            ret_list = [[7],[12],[15], [20], [25], [30], [35], [40], [45], [50]]
            #ret_list = [[10]]
        if num_dim == 2:
            ret_list = [[10, 3], [20, 3], [30, 3], [40, 3],
                        [50, 3], [50, 4],[50, 5], [50, 6]]
            #et_list = [[50, 5]]
        if num_dim == 3:
            ret_list = [[20, 3, 3], [40, 3, 3], [50, 3, 3],
                        [50, 3, 3], [50, 4, 3], [50, 5, 3], [50, 6, 3],
                        [50, 6, 3], [50, 6, 4],
                        ]
            #ret_list = [[20, 3, 3], [50, 3, 3], [50, 4, 4]]

    if unit=='per':
        """value is % from available training, not arriving at 100 means not using all the training"""
        if IS_size=='MEDIUM':
            if num_dim == 1:
                ret_list = [[10],[15], [20], [25], [30], [40], [60], [100]]
            if num_dim == 2:
                ret_list = [[10, 80],[20, 80],[40, 80],[60, 80], [80, 80],[85, 100],
                            [100, 100]]
            if num_dim == 3:
                ret_list = [[10, 80, 80], [20, 80, 80], [30, 80, 80],[50, 80, 80],[80, 80, 80]
                            ,[85, 80, 90],[90, 90, 90],[100, 100, 100]
                            ]

        if IS_size=='BIG':
            if num_dim == 1:
                ret_list = [[10],[14],[18], [20], [25], [30], [40], [60],[80], [100]]

            if num_dim == 2:
                ret_list = [[10, 25],[15, 35],[20, 45],[30, 60], [40, 70],[50, 80], [60, 90],[80, 100],
                            [100, 100]]
            if num_dim == 3:
                ret_list = [[10, 60, 60],[15, 70, 60],[20, 70, 70], #lestt than 50 in 2nd and 3rd dim will produce a knot vector error
                            [30, 80, 70], [40, 80, 80], [50, 90, 80],[65, 90, 90],[80, 100, 100],[100, 100, 100]
                            ]

    return {unit:ret_list, 'serma_supp':'Y'}

def user_app_input(num_dim,study=None):
    """Typical user data for support budget
    options:
    'fg_full_bsp_O3_splop_',
    """
    if study!=None:
        apps=aux_case_helper(num_dim,study)

    return apps


def aux_case_helper(num_dim,study):
    apps=OrderedDict()
    if study=='_fgmxsrmxs_':
        pass

    if '_XEtest_'  in study:
        n='_sct_mxsal_ksp1_'
        apps[n]={}
        apps[n]['types']=['_km_C_' for _ in range(num_dim)]
        apps[n]['k']=['_sp_' for _ in range(num_dim)]
        apps[n]['supp_type']='_sct_mxsal_'
        apps[n]['degree']=[1 for _ in range(num_dim)]
        n='_sct_rmxsal_ksp1_'
        apps[n]={}
        apps[n]['types']=['_km_C_' for _ in range(num_dim)]
        apps[n]['k']=['_sp_' for _ in range(num_dim)]
        apps[n]['supp_type']='_sct_rmxsal_'
        apps[n]['degree']=[1 for _ in range(num_dim)]

    if '_BestKM_'  in study:
        n='fg_bsp_linear_'
        apps[n]={}
        apps[n]['types']=['_bsp_' for _ in range(num_dim)]
        apps[n]['orders']=[2 for _ in range(num_dim)]
        apps[n]['supp_type']='_fg_'
        apps[n]['t']=['_linear_' for _ in range(num_dim)]

        n='fg_bsp_preplop_3_'
        apps[n]={}
        apps[n]['types']=['_bsp_' for _ in range(num_dim)]
        apps[n]['orders']=[3 for _ in range(num_dim)]
        apps[n]['supp_type']='_fg_'
        apps[n]['t']=['_preplop_' for _ in range(num_dim)]

        n='fg_bsp_oddend_3_'
        apps[n]={}
        apps[n]['types']=['_bsp_' for _ in range(num_dim)]
        apps[n]['orders']=[3 for _ in range(num_dim)]
        apps[n]['supp_type']='_fg_'
        apps[n]['t']=['_oddend_' for _ in range(num_dim)]

        n='_sct_rixsalTOL_ksp1_'
        apps[n]={}
        apps[n]['types']=['_km_C_' for _ in range(num_dim)]
        apps[n]['k']=['_sp_' for _ in range(num_dim)]
        apps[n]['supp_type']='_sct_rixsalTOL_'
        apps[n]['degree']=[1 for _ in range(num_dim)]

        n='_sct_rixsalTOL_ksp2_'
        apps[n]={}
        apps[n]['types']=['_km_C_' for _ in range(num_dim)]
        apps[n]['k']=['_sp_' for _ in range(num_dim)]
        apps[n]['supp_type']='_sct_rixsalTOL_'
        apps[n]['degree']=[2 for _ in range(num_dim)]

        n='_sct_rixsalTOL_ksp2_regul11_'
        apps[n]={}
        apps[n]['types']=['_km_C_' for _ in range(num_dim)]
        apps[n]['k']=['_sp_' for _ in range(num_dim)]
        apps[n]['supp_type']='_sct_rixsalTOL_'
        apps[n]['degree']=[2 for _ in range(num_dim)]
        apps[n]['rcond']=[1E-11 for _ in range(num_dim)]

        n='_sct_rixsalTOL_ksp2_pcall_'
        apps[n]={}
        apps[n]['types']=['_km_C_' for _ in range(num_dim)]
        apps[n]['k']=['_sp_' for _ in range(num_dim)]
        apps[n]['supp_type']='_sct_rixsalTOL_'
        apps[n]['degree']=[2 for _ in range(num_dim)]
        apps[n]['precd']='_sqrt(all)_'

        n='_sct_rixsalTOL_ksp2_pcall_regul11_'
        apps[n]={}
        apps[n]['types']=['_km_C_' for _ in range(num_dim)]
        apps[n]['k']=['_sp_' for _ in range(num_dim)]
        apps[n]['supp_type']='_sct_rixsalTOL_'
        apps[n]['degree']=[2 for _ in range(num_dim)]
        apps[n]['precd']='_sqrt(all)_'
        apps[n]['rcond']=[1E-11 for _ in range(num_dim)]

    if len(apps.keys())<1:
        raise ValueError('no apps given by the user')
    if len(apps.keys())!=len(set(apps.keys())):
        raise ValueError('apps repetition')
    return apps


def closest_serma_disct(dim, vec, less_sermapoints=None):
    """
    generating thetandSERMA support, baon hdoisitions.
    """

  # Hardcoded "andavalues for a UO2
    serma_dic = {}
    serma_dic['BURNUP'] = [0, 9.375, 18.75,  37.5,   75,   150,237.5,  325,  412.5,   500,   625,   750,  1000,
                               1250,  1500,  1750,  2000, 2500, 3000, 3500,   4000,  4500,  5000,  5500,  6000,
                               6500,  7000,  7500,  8000, 8500, 9000, 10000, 12000, 15000, 18000, 20000, 23000,
                              26000, 30000, 35000, 40000, 4.48734492e+04]
    serma_dic['TEMPERATURE_COMBUSTIBLE'] = [295, 303, 650, 800, 1200, 2000]
    serma_dic['CONCENTRATION_EN_BORE'] = [0, 300, 600, 900, 1600]

    # print dim
    if dim.dim_name not in serma_dic:
        raise RuntimeError("The name of the %s is cannot be found un hard-coded %s" %
                           (dim.dim_name, serma_dic.keys()))

    # Swhich this off with a denser grid!
    if dim.dim_name == 'BURNUP' and less_sermapoints == True:
        serma_dic['BURNUP'].pop(2)
        serma_dic['BURNUP'].pop(10)
        serma_dic['BURNUP'].pop(11)
        serma_dic['BURNUP'].pop(17)
        serma_dic['BURNUP'].pop(20)
        serma_dic['BURNUP'].pop(24)
        serma_dic['BURNUP'].pop(30)
    # print len(serma_dic['BURNUP'])

    serma_supp = []
    serma_supp_index = []
    for serma_point in serma_dic[dim.dim_name]:
        # for each serma_point the whole domain is run to find the difference but
        # this is not required per se
        diff = [abs(serma_point - val) for val in vec]
        i = np.argwhere(diff == min(diff))[0][0]
        serma_supp.append(vec[i])
        serma_supp_index.append(i)

    return sorted(set(serma_supp_index))


