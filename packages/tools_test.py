# !/usr/bin/env python
"""
This file containes some important tests.
"""

from packages.tools_SG_common import *

def TT_separation_test(TT_dic):
    for supp in TT_dic['supp']:
        # print TT_dic['supp'][supp].__doc__
        checking_TT_separation(
            list(
                TT_dic['supp'][supp].grid_nmlz_idx), list(
                TT_dic['test']['user_level'].grid_nmlz_idx))
    for supp in TT_dic['supp']:
        checking_TT_separation(list(TT_dic['supp'][supp].xs('_1_U235_abso_2')), list(TT_dic[
            'test']['user_level'].xs('_1_U235_abso_2')))

def true_value_test(app_dic,TT_dic,cont_UO2medium_1,complete_specialized_iso=False,verbose=False):
    """
    The objective of this test is to analyze if true values are met for k_inf when evaluating in the support for an interpolation scheme
    These properites are tested and do not belong to standard output as some require estimation of concentration and this is not a standard feature.
    """
    for app_n,app_v in app_dic.iteritems():
        #app_v=app_dic['f_bsp_linear_dg1_fg_n51']
        if 'f_bsp_linear_dg1' in app_n:
            do_true_val_test(app_v,TT_dic,cont_UO2medium_1,complete_specialized_iso,verbose)


def do_true_val_test(app_v,TT_dic,cont_UO2medium_1,complete_specialized_iso,verbose):
    app_v.variables_mungling1()
    app_v.variables_mungling2() 
    if verbose==True: print 'for ', app_v

    # finding isotopes of interest for every xs
    xs_interest={}
    for mxs_name, mxs_tup in app_v.get_prp('a2','mxs_nametuple_space').iteritems():
        xs_interest[mxs_name]=app_v.find_iso_of_mxs(mxs_tup)

    if complete_specialized_iso==True:
        for mxs_n,xs_list in xs_interest.iteritems():
            if 'nufi' not in mxs_n:
                if len(xs_list)!=len([i for i in cont_UO2medium_1.XIRG['I']['name'] if i!='MACRT']):
                    raise ValueError('isotopes not properly accounted for')

    #analyzing mxs error in true value points
    supp_grid=list(app_v.get_prp('supp','grid_nmlz_val'))
    mxs_in_supp={}
    mxs_in_supp['true']={}
    mxs_in_supp['estimated']={}
    for mxs_n,xs_list in xs_interest.iteritems():
        xs_sigma_conc_dict={}
        xs_sigma_conc_dict['sigma_true']={}
        xs_sigma_conc_dict['sigma_estimated']={}
        xs_sigma_conc_dict['conc']={}
        for xs in xs_list:
            xs_sigma_conc_dict['sigma_true'][xs]=list(app_v.get_fprp('supp','xs',xs)) # This is the 'true value'
            xs_sigma_conc_dict['conc'][xs]=list(app_v.get_fprp('supp','conc_grid',xs))
            _,f=app_v.get_xs_app(xs)
            xs_sigma_conc_dict['sigma_estimated'][xs]=[f(x) for x in supp_grid]
        mxs_in_supp['true'][mxs_n]=app_v.mxs_calcul(xs_sigma_conc_dict['sigma_true'],xs_sigma_conc_dict['conc'])
        mxs_in_supp['estimated'][mxs_n]=app_v.mxs_calcul(xs_sigma_conc_dict['sigma_estimated'],xs_sigma_conc_dict['conc'])
        if verbose==True: print ".. comparing MXS build using sigma values and sigma_ch for support points"
        if max([abs(s-sch) for s,sch in itertools.izip_longest(mxs_in_supp['true'][mxs_n], mxs_in_supp['estimated'][mxs_n], fillvalue=None)])>1E-16:
            raise ValueError('mxs error in support is too high')

    # comparing to MACRT ? not present in the majority of standarized data sets
    if 'MACRT' in cont_UO2medium_1.XIRG['I']['name']:
        mxs_in_supp['MACRT']={}
        for mxs in mxs_in_supp['true']:
            aux_name_sale='_1_MACRT_'+mxs.split('_')[2]+'_'+mxs.split('_')[3]
            mxs_in_supp['MACRT'][mxs]=list(app_v.get_fprp('supp','xs',aux_name_sale))
            mxs_rebuilt=mxs_in_supp['true'][mxs]
            if verbose==True: print ".. comparing rebuild MXS in supp point with actaul MACRT values"
            error_aux=max([abs(s-sch) for s,sch in itertools.izip_longest(mxs_in_supp['MACRT'][mxs], mxs_rebuilt, fillvalue=None)])*100/max(mxs_in_supp['MACRT'][mxs])
            if error_aux>1E-4:
                print mxs, 'failed test of rebuild. This is normal for nufi 1 and 2 if the negative values of nufi RES are not being considered...'
                print 'in this case RES.keys()=',cont_UO2medium_1._a2_data['xs']['1']['MACR'].keys(), error_aux


    required_vec=['_1_abso_1', '_1_tran012_1', '_1_nufi_1', '_1_tran021_2', '_1_nufi_2', '_1_abso_2']
    if all([req in xs_interest.keys() for req in required_vec]):

        k_true=app_v.k_calcul('k_inf_classic',mxs_in_supp['true'])
        k_esti=app_v.k_calcul('k_inf_classic',mxs_in_supp['estimated'])

        # analyzing that interpolation actually preserves k values
        for key in k_true:
            if verbose==True: print '.. comparing k between rebuild MXS using true value and estimated interpolation values in support'
            k_interpol_error=max([1E5*abs(k-k_e) for k,k_e in itertools.izip_longest(k_true[key], k_esti[key], fillvalue=None)])
            if k_interpol_error>1:
                raise ValueError('k calculated in true value has too high of an error')

        # analyzing that at least for the 'rebuild' k_inf support point reproduce true values
        if verbose==True: print '.. comparing k between rebuild MXS using estimated interpolation values in support and rebuild k in a2 (using true micro xs)'
        k_rebuild_error=max([1E5*abs(k-k_e) for k,k_e in itertools.izip_longest(app_v.get_fprp('supp','k_inf','_1_k_inf_rb'), k_esti['_1_k_inf_rb'], fillvalue=None)])
        if k_rebuild_error>1:
            raise ValueError('estimated k from reconstructed xs failes to reproduce reconstructed k inf in a2')

        # analyzing a2 k inf using MACRT xs
        if 'MACRT' in cont_UO2medium_1.XIRG['I']['name']:
            k_macrt=app_v.k_calcul('k_inf_classic',mxs_in_supp['MACRT'])
            if verbose==True: print '.. comparing k between rebuild k values using build MXS in true point and MACRT'
            k_rebuild_error=max([1E5*abs(k-k_e) for k,k_e in itertools.izip_longest(k_esti['_1_k_inf_rb'], k_macrt['_1_k_inf_rb'], fillvalue=None)])
            if k_rebuild_error>1:
                print 'WARNING k estimated in interpolation with MXS rebuild and k estimated with MACRT differ',k_rebuild_error
                print 'this can be corrected by including negative nufi res values'
            if verbose==True: print '.. comparing k between a2 values  and rebuild using MACRT'
            k_real_error=max([1E5*abs(k-k_e) for k,k_e in itertools.izip_longest(app_v.get_fprp('supp','k_inf','_1_k_inf_rb'), k_macrt['_1_k_inf_rb'], fillvalue=None)])
            if k_real_error>1:
                print 'WARNING: k inf obtained from a2 file differs from the obtained with the formula using MACRT. in pcm:',k_real_error
        if verbose==True: print '.. comparing k between a2 values  and rebuild using MACRT'
        k_real_error=max([1E5*abs(k-k_e) for k,k_e in itertools.izip_longest(app_v.get_fprp('supp','k_inf','_1_k_inf_a2'), k_esti['_1_k_inf_rb'], fillvalue=None)])
        print 'rebuild vs [pcm] a2=',k_real_error
        if k_real_error>1:
            print 'WARNING: k inf obtained from a2 file differs from the obtained with reconstruction with estimated k in supp points. in pcm:',k_real_error
        app_v.destroy_aux_f()



def checking_TT_separation(supp_data, test_data):
    """
     Cheking that intersection of support and test is empty
    """

    if isinstance(supp_data[0], type(tuple)):
        for supp_point in supp_data:
            if supp_point in test_data:
                print supp_point
                raise ValueError('supp point in test')
        return 0
    if isinstance(supp_data[0], type(np.float64)):  # !# is this ok?
        sign_cifr = 9
        # print round(supp_point, sign_cifr)
        aux = [round(g, sign_cifr) for g in test_data]
        for supp_point in supp_data:
            if round(supp_point, sign_cifr) in aux:
                print supp_point
                raise ValueError('supp point in test')
        return 0
        # sys.exit()

def supp_septaration_test(app_dic,TT_dic):
    """
    Check proper indexation of grid points.
    """
    test_idx=list(TT_dic['test']['user_level'].get_prp_valid('grid_nmlz_idx'))
    test_val=list(TT_dic['test']['user_level'].get_prp_valid('grid_nmlz_val'))
    idx_breach=[]
    val_breach=[]
    for app in app_dic.values():
        #print app
        for point in app.get_prp('supp','grid_nmlz_idx'):
            idx_breach.append(point in test_idx)
        for point in app.get_prp('supp','grid_nmlz_val'):
            val_breach.append(point in test_val)
    if idx_breach!=val_breach:
        raise RuntimeError('Mayor problem, lost coherence between indexation and values')
    if True in idx_breach:
         raise RuntimeError('support points in test detected!')
