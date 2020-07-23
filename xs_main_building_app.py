# --*-- coding:utf-8 --*--
"""
This file generates the approximations for the considered supports.
* defines supports and TT dicts
* Defines approximations using the clases
* Manages the supervised learning if any
"""

from __future__ import division
from packages.data_clases import *
from packages.tools_test import *
from packages.tools_SG_common import *
from packages.tools_IO import *
import packages.supervised_learning as spv
from IPython import embed

def build_support(num_dim,cont_UO2medium_1,supports):
    ret_dict=OrderedDict()
    for supp in supports:
        n=str(np.prod([len(t)for t in supp]))
        #print ' ... building support n='+str(n)
        name = 'fg_' + 'd' + str(num_dim) + '_n' + str(n) + '_'
        A = Cartesian_support(supp, cont_UO2medium_1, name)
        A.Cartesian_2_iterator_core()
        ret_dict[name] = A
    return ret_dict

def master_initialization(usedic):
	"""
	General data init
	"""
	num_dim=usedic['num_dim']
	pth =usedic['app']['root_I']+usedic['A2_case']+usedic['app']['f_I']
	name = usedic['run_case']
	apollo2_data, a2_meta = load_obj_pickle(pth, name,load_type='serial')
	a2_meta.path_app_O =usedic['app']['root_I']+usedic['A2_case']+usedic['rlt']['f_I']
	a2_meta.app_O_n = usedic['run_case']+usedic['study']
	cont_UO2medium_1 = ApolloContainer(apollo2_data, num_dim)
	print '... computing mxs' #required as only xs and concentration are present in a2 data
	cont_UO2medium_1.compute_mxs()
	cont_UO2medium_1.compute_k_inf_rb('k_inf_classic') #k_inf from a2 may be somewhat trubelsome
	I=cont_UO2medium_1.xs_importance_wrt_mxs(key='at_grid')
	I_flat={xs_n:xs_v for mxs in I.keys() for xs_n,xs_v in I[mxs].items()}
	#print I_flat.keys()
	a=cont_UO2medium_1.mxs_nametuple_space.values()[0]
	AL_LIM=0.2
	valid_apriori_xs_namespace=[]
	for mxs_tuple in cont_UO2medium_1.mxs_nametuple_space.values():
	    xs_by_imp=[]
	    for xs in cont_UO2medium_1.find_iso_of_mxs(mxs_tuple):
	        if 'MACRT' in xs: continue
	        #xs_by_imp.append([xs,np.average(I_flat[xs])])
	        xs_by_imp.append([xs,max(I_flat[xs])])
	    #xs_by_imp=list(reversed(sorted(xs_by_imp,key=lambda x:x[1])))
	    valid_apriori_xs_namespace.append([name_integratedImp[0] for name_integratedImp in xs_by_imp if name_integratedImp[1]>AL_LIM])
	valid_apriori_xs_namespace=np.concatenate(valid_apriori_xs_namespace)
	# embed()
	print '... computing training'
	cont_UO2medium_1.compute_train_gen_bybudget(
	     training_budget=get_user_TT_division(usedic['IS_size'],num_dim), less_sermapoints=usedic['less_sermapoints'])
	print 'total training points=',len(cont_UO2medium_1.train_grid_idx)
	print ' ... building test'
	TT_dic = OrderedDict()
	TT_dic['test'] = OrderedDict()
	A = Test(cont_UO2medium_1, cartesian_training_index=cont_UO2medium_1.train_cartesian_idx)
	A.Cartesian_2_iterator_core()
	if '_testtriche_' in usedic['study']:
	    flag_triche='YES'
	else:
	    flag_triche='NO'
	A.sparsefy_iterator_core(usedic['density_regulizer'],test_lim=10000,bound_to_training='yes',triche_grid=flag_triche)
	test_grid=list(A.get_prp_valid('grid_nmlz_val'))
	if len(test_grid)<2:
	    raise ValueError('U run out of test points')
	TT_dic['test']['user_level'] = A
	# Building support generator
	supports_budget = get_user_supp_input(usedic['IS_size'],num_dim)
	print '... building support for '
	cont_UO2medium_1.compute_supp_gen(cp.deepcopy(supports_budget))
	supports = cont_UO2medium_1.supp_gen_idx
	TT_dic['supp'] = build_support(num_dim,cont_UO2medium_1,supports)
	supp_resources={}
	supp_resources['I']=I
	supp_resources['I_flat']=I_flat
	supp_resources['valid_apriori_xs_namespace']=valid_apriori_xs_namespace
	supp_resources['supports']=supports
	return num_dim, apollo2_data,a2_meta,cont_UO2medium_1,supp_resources,TT_dic,usedic

def build_app(cont_UO2medium_1,supp_v,app_n,app_v,test_grid,I_flat,valid_apriori_xs_namespace,previous_idx=None,verbose=False):
    """
    note that supp_v is cartesian. If FG this is retained, otherwise the equivalent bufget is obtained and used to define a sacterred support
    """
    if '_fg_'in app_v['supp_type']:
        supp_v=supp_v
        log_al_app={}
    if '_sct_'in app_v['supp_type']:
        budget=supp_v.budget
        #sys.exit()
        supp_v,log_al_app=spv.supp_scattered(cont_UO2medium_1,I_flat,budget,app_v,valid_apriori_xs_namespace,test_grid_idx=test_grid,previous_idx=previous_idx,verbose=verbose)
    #embed()
    A=spv.instanziate_app(app_n,app_v,supp_v)
    #A.approximate_xs() #repeated!
    A._test_inst=test_grid
    return A,A.dict_special_app['app_name'],log_al_app


def serial_main(usedic):
    """
    Standarization to convenient a2-python structure
    """

    def build_app_dict():

        apps = user_app_input(num_dim,study=usedic['study'])
        apps_loop={}
        apps_nloop={}
        for app_n,app_v in user_app_input(num_dim,study=usedic['study']).iteritems():
            if 'wal_' not in app_v['supp_type']:
                apps_loop[app_n]=app_v
            else:
                apps_nloop[app_n]=app_v
        app_dic = OrderedDict()
        log_al=OrderedDict()
        verbose=False
        for supp_v in TT_dic['supp'].values():
            for app_n,app_v in apps_loop.iteritems():
                previous_idx=find_previous_app(app_v,app_dic)
                print '... building xs approximation for', supp_v.name, app_n
                A,app_save,log_al_app=build_app(cont_UO2medium_1,supp_v,app_n,app_v,TT_dic['test']['user_level'],supp_resources['I_flat'],supp_resources['valid_apriori_xs_namespace'],previous_idx=previous_idx,verbose=verbose)
                try:
                    log_al[app_save]=log_al_app
                except UnboundLocalError:
                    pass
                app_dic[app_save] = A
                #print app_dic.keys()

        #sys.exit()
        if apps_nloop:
            #print 'here'
            for app_n,app_v in apps_nloop.iteritems():
                print '... building xs approximation for', app_n
                A,app_save,log_al_app=build_app(cont_UO2medium_1,supp_v,app_n,app_v,TT_dic['test']['user_level'],supp_resources['I_flat'],supp_resources['valid_apriori_xs_namespace'],verbose=verbose)
                try:
                    log_al[app_save]=log_al_app
                except UnboundLocalError:
                    pass
                app_dic[app_save] = A
        #sys.exit()
        return app_dic, log_al
    # action of standarization
    beg=time.time()
    num_dim, apollo2_data,a2_meta,cont_UO2medium_1,supp_resources,TT_dic,usedic=master_initialization(usedic)
    app_dic,log_al = build_app_dict()
    print 'a total jobs=',len(app_dic.keys()),'was run in processes N=1', 'in a time of [min]',(time.time()-beg)/60
    save_obj_pickle([app_dic,log_al,cont_UO2medium_1,supp_resources['I_flat'],a2_meta,TT_dic], a2_meta.path_app_O, a2_meta.app_O_n)
    print 'sucess'
    print 'density regul=', usedic['density_regulizer']


# permits to lunch runs of calculation
if __name__ == "__main__":
    #print " ... building approximations"
    IS_size_run=['BIG']
    dim_run=[1]
    OS_run=['_X1I31R6G2_']
    study_run=['_BestKM_']
    run_case=[]
    for IS_size in IS_size_run:
        for num_dim in dim_run:
            for OS_size in OS_run:
                for study in study_run:
                    aux=user_case()
                    aux['num_dim']=num_dim
                    aux['OS_size']=OS_size
                    aux['study']=study
                    aux['IS_size'] = IS_size
                    aux['run_case'] = aux['A2_case'] +aux['OS_size'] + aux['IS_size']+'_' + 'd' + str(aux['num_dim'])
                    print aux['run_case']
                    run_case.append(aux)
    print len(run_case)
    #run_case=run_case[sys.argv[1]] # to run from user case
    #run_case=[user_case()]
    goto=0
    #goto=6
    for case in run_case:
    	random.seed(a=0)
    	serial_main(case)
