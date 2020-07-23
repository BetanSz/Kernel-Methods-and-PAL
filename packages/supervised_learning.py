# --*-- coding:utf-8 --*--
"""
This file containes active learning functions
"""

from packages.tools_SG_common import *
from packages.data_clases import *
from IPython import embed

def instanziate_app(app_n,app_v,sup_init,acc=False,user_xs_names=False):
    A=None
    if all(['_bsp_' in type_i for type_i in app_v['types']]):
        A = Approximation_SPL(app_n,app_v,sup_init)
    if all(['_pp_' in type_i for type_i in app_v['types']]):
        A = Approximation_SPL(app_n,app_v,sup_init)
    if all(['_km_C_' in type_i for type_i in app_v['types']]):
        A = Approximation_CKM(app_n,app_v,sup_init)
    if all(['_km_sk_' in type_i for type_i in app_v['types']]):
        A = Approximation_SK(app_n,app_v, sup_init)
    if not A:
        raise ValueError('cant find approximation')

    A.approximate_xs(user_xs_names=user_xs_names)
    if acc!=False:
        A.accelerate(acc)
    return A


def supp_scattered(data_set,I_flat, budget,app_v,valid_apriori_xs_namespace,test_grid_idx=None,previous_idx=None, verbose=False):
    """
    Allocates a scatter supp
    Each AL process is independent of the previous one. Though using older points would accelerate allot the whole thing
    """
    if verbose==True: print 'Retriving data_set..'
    dict_al=OrderedDict()
    def unloading_shell():
        return data_set.train_grid_idx,data_set.grid_nmlz_val
    training_idx,A2_X=unloading_shell()

    def log_learning(dict_supp_sct,dict_al,N):
        dict_supp_sct['novo_point']=[data_set.grid_nmlz_val[val] for val in dict_supp_sct['novo_idx']]
        dict_al[str(N)]=dict_supp_sct
    if '_rdm_' in app_v['supp_type']:
        N_init=budget
    if 'al' in app_v['supp_type']:
        N_init=3

    if verbose==True: print '...generating random supp'
    dict_supp_sct=supp_random_learning(training_idx, N_init)
    if previous_idx:
        supp_idx=previous_idx
    else:
        supp_idx=dict_supp_sct['novo_idx']

    sup_v=Scattered_support(data_set,'random',valid_idx=supp_idx)
    log_learning(dict_supp_sct,dict_al,len(supp_idx))
    if verbose==True: print 'starting adding points...,',len(supp_idx),budget,N_init,app_v
    if 'wal_' not in app_v['supp_type']:
        while len(supp_idx)<budget:
            beg=time.time()
            dict_supp_sct=supp_active_learning(data_set,I_flat,training_idx,sup_v,app_v,supp_idx,valid_apriori_xs_namespace,safe_flag=False)
            if verbose==True: print 'doing supp AL',time.time()-beg,'for ',len(supp_idx)-budget
            beg=time.time()
            supp_idx=sorted(supp_idx+dict_supp_sct['novo_idx'])
            sup_v=Scattered_support(data_set,'iteration in al process',valid_idx=supp_idx)
            log_learning(dict_supp_sct,dict_al,len(supp_idx))
    else:
        dict_xsw_supp_sct=xsw_supp_active_learning(data_set,I_flat,training_idx,sup_v,app_v,verbose=verbose)
        sup_v=Scattered_support(data_set,'iteration in al process',valid_idx=dict_xsw_supp_sct,)
        dict_al={}

    if verbose==True: print '...finish adding points'
    if len(supp_idx)!=len(set(supp_idx)):
        raise ValueError('AL reapeting points')
    return sup_v,dict_al

def xsw_supp_active_learning(A2_data,I_flat,training_idx,sup_init,app_v,verbose=False):
    """
    Here the termination criterion is to achive an error smaller than TOL for every XS separadly
    Then it is also consider that error but wheighted with the improtance allowing big error of uniportant XSs
    """
    def unified_idx_valid(data,training_idx,excluded_idx):
        '''
        Returns an index list always restpect to the training_idx
        '''
        aux=[[idx,val] for idx,val in enumerate(data) if idx in training_idx and idx not in excluded_idx]
        return zip(*aux)

    # print app_v
    # sys.exit()
    a2_data=A2_data.grid_nmlz_val
    TOL=5E-2
    if 'TOL' in app_v:
        TOL=app_v['TOL']
    xsw_supp={}
    #verbose=True
    for xs_n in A2_data.xs_namespace:
        supp_idx=cp.deepcopy(sup_init.valid_idx)
        Ltraining=len(unified_idx_valid(a2_data,training_idx,supp_idx)[0])
        Lsupp=len(supp_idx)
        error=1E10
        y=A2_data.xs(xs_n)
        #if
        if app_v['supp_type']=='_sct_xsiwal_':
            modifier=np.mean(I_flat[xs_n])
        if app_v['supp_type']=='_sct_xswal_':
            modifier=1.0
        if verbose==True: print xs_n,modifier
        while error>TOL and Lsupp<Ltraining:
            supp_v=Scattered_support(A2_data,'iteration in al process',valid_idx=supp_idx)
            x_idx_vec,x_valid=unified_idx_valid(a2_data,training_idx,supp_idx)
            acc=False
            A=instanziate_app('iteration',app_v,supp_v,acc=acc,user_xs_names=[xs_n]) # at least I only work in the XS of interest, even if there is some recalculation of K_training
            _,f=A.get_xs_app(xs_n)
            xs_in_training=[f(x_v,idx=x_i) for x_v,x_i in itertools.izip_longest(x_valid,x_idx_vec,fillvalue=None)]
            _,y_in_training=unified_idx_valid(y,training_idx,supp_idx)
            errors=[[x_idx,abs(y_ch/yi-1)*100] for x_idx,y_ch,yi in itertools.izip_longest(x_idx_vec,xs_in_training,y_in_training, fillvalue=None)]
            error=np.mean([err[1] for err in errors])* modifier
            idx_novo=max(errors, key=lambda x:x[1])[0]
            supp_idx.append(idx_novo)
            supp_idx=sorted(supp_idx)
            Lsupp=len(supp_idx)
        if verbose==True: print len(supp_idx), len(training_idx)
        xsw_supp[xs_n]=supp_idx

    return xsw_supp

def supp_active_learning(A2_data,I_flat,training_idx,sup_init,app_v,supp_idx,valid_apriori_xs_namespace,safe_flag=True):
    """
    Performes AL.
    index novo is a list, as in the future batch AL may ve available
    """
    training_eff=[val for val in training_idx if val not in supp_idx] #if valid not in supp_idx (but in training and not in test)
    def valid(data):
        """
        This must bu run in and index wrt the total A2 data
        """
        return [data[val] for val in training_eff]
    x=A2_data.grid_nmlz_val

    if safe_flag:
        if len(valid(x))!=len(training_idx)-len(supp_idx):
            raise ValueError('Effective training not well defined')
    acc=False
    A=instanziate_app('iteration',app_v,sup_init,acc=acc)
    # Defining the xs namespace uses to label supp points
    namespace=A.get_prp('a2','xs_namespace')
    if '_Uxsal_' in app_v['supp_type']:
        xs_n='_1_U235_nufi_2'
        namespace=[xs_n]
        I_flat_copy={}
        I_flat_copy[xs_n]=[1 for _ in x]
    if 'TOL' in app_v['supp_type']:
        namespace=valid_apriori_xs_namespace

    def lossF_calcul():
        """
        This functions gives back the dict of errors
        """
        #print valid(x)
        x_idx,x_valid=zip(*list(enumerate(valid(x))))
        # Building xs
        xs_in_training={}
        for xs_n in namespace:
            if 'MACRT' in xs_n:continue
            if 'MACR' in xs_n and 'nufi' in xs_n: raise ValueError('you may want to check what XS are being consider in AL')
            _,f=A.get_xs_app(xs_n)
            xs_in_training[xs_n]=[f(x_v,idx=x_i) for x_v,x_i in itertools.izip_longest(x_valid,x_idx,fillvalue=None)]

        errors={}
        if 'xsal'in app_v['supp_type']:
            if '_ixsal_' in app_v['supp_type'] or '_Uxsal_' in app_v['supp_type'] or '_misalTOL_' in app_v['supp_type']:
                error_type=lambda y_ch,y,I: abs(y_ch-y)*I
            if '_rixsal_' in app_v['supp_type'] or '_rUxsal_' in app_v['supp_type'] or '_rixsalTOL_' in app_v['supp_type']:
                error_type=lambda y_ch,y,I: abs(y_ch/y-1)*I
            if '_rnixsal_' in app_v['supp_type']:
                error_type=lambda y_ch,y,Ii: abs(y_ch/y-1)
            for xs_n in xs_in_training.keys():
                errors[xs_n]=[error_type(y_ch,yi,Ii) for y_ch,yi,Ii in itertools.izip_longest(xs_in_training[xs_n],valid(A2_data.xs(xs_n)),valid(I_flat[xs_n]),fillvalue=None)]

        if 'M' in app_v['supp_type'] or 'k' in app_v['supp_type']:
            error_type=lambda y_ch,y: abs(y_ch/y-1)
            A.variables_mungling1() # allows to use find_iso_of_mxs generating the self. name spaces
            xs_interest={}
            for mxs_n, mxs_tup in A.get_prp('a2','mxs_nametuple_space').iteritems():
                xs_interest[mxs_n]=A.find_iso_of_mxs(mxs_tup) #what xs are required for that mxs
            # Calculating mxs
            mxs_in_training={}
            for mxs_n in xs_interest:
                conc_dict=OrderedDict()
                xs_dict=OrderedDict()
                for xs_n in xs_interest[mxs_n]:
                    conc_dict[xs_n]=valid(A2_data.conc_grid(xs_n))
                    xs_dict[xs_n]=xs_in_training[xs_n]
                mxs_in_training[mxs_n]=A.mxs_calcul(xs_dict,conc_dict)

        if '_Mal_' in app_v['supp_type']:
            errors={}
            for mxs_n,mxs_v in mxs_in_training.iteritems():
                errors[mxs_n]=[error_type(y_ch,yi) for y_ch,yi in itertools.izip_longest(mxs_v,valid(A2_data.mxs(mxs_n)),fillvalue=None)]

        if '_kal_' in app_v['supp_type']:
            required_vec=['_1_abso_1', '_1_tran012_1', '_1_nufi_1', '_1_tran021_2', '_1_nufi_2', '_1_abso_2']
            if all([req in xs_interest.keys() for req in required_vec]):
                k_si=A.k_calcul('k_inf_classic',mxs_in_training)
                k_training=k_si['_1_k_inf_rb']
                errors['k']=[error_type(y_ch,yi) for y_ch,yi in itertools.izip_longest(k_training,valid(A2_data.k_inf('_1_k_inf_rb')),fillvalue=None)]
            else:
                raise ValueError('AL on k but MXS required not present')

        return errors
    errors=lossF_calcul()
    def max_xs_calcul():
        return [[xs_n,errors[xs_n].index(max(errors[xs_n])),max(errors[xs_n])] for xs_n in errors]
    def selected_xs_calcul():
        return max(max_xs,key=lambda x:x[2])
    max_xs=max_xs_calcul()
    selected_xs=selected_xs_calcul()
    position_max=selected_xs[1]
    return {'novo_idx':[training_eff[position_max]],'max_xs':max_xs,'selected_xs':selected_xs,'type':app_v['supp_type']}

def supp_random_learning(training_idx, budget):
    """
    the only relevant part here is
    random.choice(training_idx)
    this can be fusioned with normal al
    """
    supp_idx=[]
    while len(supp_idx)<budget:
        supp_idx.append(random.choice(training_idx))
        supp_idx=list(set(supp_idx))
        supp_idx=sorted(supp_idx)
    return {'novo_idx':supp_idx,'type':'random'}
