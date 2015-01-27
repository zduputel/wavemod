'''
My Utils

Written by Z. Duputel, September 2013
'''

# Externals
from os.path import exists,isdir
from numpy   import ndarray
import shutil
import os  

def rm(ifiles):
    '''
    Remove file(s) or directory(ies)
    Args:
         ifiles:  file, directory or list of files, directories
    '''

    if type(ifiles)!=list and type(ifiles)!=ndarray:
        ifiles = [ifiles]
    for ifile in ifiles:
        if exists(ifile):
            if isdir(ifile):
                shutil.rmtree(ifile)
            else:
                os.remove(ifile)
                
    # All done
    return


def mkdir(idir,pflag=False):
    '''
    Create directory
    Args:
         idir:   directory to be created
         pflag: (optional, default=False)
               False: will     remove directory if exists
               True:  will NOT remove directory if exists
    '''        

    if pflag and not exists(idir):
        os.mkdir(idir)
    else:
        rm(idir)
        os.mkdir(idir)
        
    # All done
    return


def parse_config(cfg_file):
    '''
    Parse my config files
    Args:
         cfg_file: configuration filename
    '''

    config = {}
    assert exists(cfg_file), 'Cannot read %s: No such file'%(cfg_file)
    try:
        config_lines = open(cfg_file, 'r').readlines()
        for line in config_lines:
            if line.find('#')==0:
                continue
            if line.rstrip():
                key,value = line.strip().split(':')
                key   = key.strip()
                value = value.strip()
                if config.has_key(key):
                    config[key].append(value)
                else:
                    config[key]=value
    except:
        stderr.write('Error: format  %s\n'%cfg_file)
        exit(1)

    # All done
    return config
