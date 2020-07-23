# --*-- coding:utf-8 --*--
"""
Functions that control IO
"""

import cPickle
from packages.tools_SG_common import *
#import dill


def worker_save(obj, adress):
    with open(adress + '.cpkl', 'wb') as f:
        cPickle.dump(obj, f, cPickle.HIGHEST_PROTOCOL)


def worker_load(adress):
    """
    sub process for loading piece of dic
    """
    with open(adress, 'rb') as f:
        return cPickle.load(f)


def valid(string, tYpe):
    """
    transforms incoming string into a valid one for saving the file or generating folders.
    """

    if tYpe == 'grl':  # general changes here
        string = string.replace('%', r'p')
        string = string.replace('1/cm', r'1_ov_cm')

    return string


def save_obj_pickle(obj, adress, name, save_type='single_file'):
    """
    Saves object as a pickle object
    """

    print 'saving pickle data...', adress + name
    grow(adress)
    if save_type == 'single_file':
        with open(adress + name + '.cpkl', 'wb') as f:
            cPickle.dump(obj, f, cPickle.HIGHEST_PROTOCOL)

    #"""
    # Serial divided
    if save_type == 'serial_divided':
        os.system("rm " + adress + "/*.cpkl" + "> /dev/null 2>&1")
        beg = time.time()
        for dic in obj:
            # separating test_dic of xs_dic
            if dic['name'] == 'xs_dic':
                with open(adress + name + '/' + dic['name'] + '.cpkl', 'wb') as f:
                    cPickle.dump(dic, f, cPickle.HIGHEST_PROTOCOL)
            else:
                # separating each test
                for test in dic['test'].keys():
                    with open(adress + name + '/' + dic['name'] + '_' + test + '.cpkl', 'wb') as f:
                        cPickle.dump(dic['test'][test], f, cPickle.HIGHEST_PROTOCOL)
        print time.time() - beg  # 0.51


def load_obj_pickle(adress, name, load_type=None, load_way=None):
    """
    Saves object as a pickle object
    """
    import cPickle

    print 'loading data...', adress + name + '.cpkl'

    if load_type == 'serial' or load_type == None:
        if load_way == None or load_way == 'cpickle':
            with open(adress + name + '.cpkl', 'rb') as f:
                return cPickle.load(f)



def saving_text(lines, results_path, file_name, subfolder=None, Format=None):
    pth = results_path
    grow(pth)
    if subfolder != None:
        pth = pth + '/' + subfolder
    grow(pth)

    if Format == None:
        extension = 'txt'
    else:
        extension = Format

    pth = pth + '/' + file_name + '.' + extension
    print 'saving ', pth
    with open(pth, 'w') as f:
        f.writelines(lines)


def grow(pth):
    """
    Unit function of generation for the tree of results
    """

    aux = []
    for folder in pth.split('/'):
        if folder == '':
            continue
        aux.append('/' + folder)
        create = ''.join(aux)
        # print create,os.path.isdir(create)
        if not os.path.isdir(create):
            print "making dir", create
            os.mkdir(create)

    if not os.path.isdir(pth):
        print "making dir", pth
        os.mkdir(pth)
        print 'here'
        sys.exit()

def save_table(results_path, head, bulk, foot, table_name, file_type):
    """
    saves the table
    """
    print 'saving table...', results_path, table_name, file_type
    grow(results_path)
    text_file = open(results_path + table_name + file_type, 'w')

    for body in [head, bulk, foot]:
        for line in body:
            print >> text_file, format(line)

    text_file.close()