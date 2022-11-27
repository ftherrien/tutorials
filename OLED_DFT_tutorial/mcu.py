import os.path
import numpy as np

def read_WAVEDER(file='WAVEDER'):
    '''Read the WAVEDER file                                                                                                                                  
                                                                                                                                                              
        the matrix of the derivative of the cell periodic part                                                                                                
        of the wavefunctions with respect to i k is:                                                                                                          
        cder = CDER_BETWEEN_STATES(m,n,1:NK,1:ISP,j)= <u_m|  -i d/dk_j | u_n> = - <u_m| r_j |u_n>                                                             
    '''
    if not os.path.isfile(file):
        print('Cannot find the %s file. Check the path:' % file)

    from scipy.io import FortranFile
    data = FortranFile(file, 'r')
    nb_tot, nbands_cder, nkpts, ispin = data.read_record(dtype= np.int32)
    nodesn_i_dielectric_function = data.read_record(dtype= np.float)
    wplasmon = data.read_record(dtype= np.float).reshape(3,3)
    cder = data.read_record(dtype= np.complex64).reshape(nb_tot, nbands_cder, nkpts, ispin, 3, order="F")
    cder = np.moveaxis(cder, [2,3],[0,1])

    return cder, nodesn_i_dielectric_function, wplasmon
