#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#


'''Script for building the cffi part manually. This scirpt is optional and can
be deleted if desired. '''

import os
import cffi

def cffibuilder(header_file, libary):
    """builds cffi libary"""

    stream = os.popen(f'cc -E {header_file}')
    output = stream.read()

    ffibuilder = cffi.FFI()
    ffibuilder.set_source("dftb_cffi",
                          f'#include "{header_file}"',
                          libraries=[f'{libary}'])

    ffibuilder.cdef(output +
                    'extern "Python" void _calc_extpot_callback(void *, double'
                    ' *, double *);\n' +
                    'extern "Python" void _calc_extpotgrad_callback(void *, '
                    'double *, double *);\n')

    ffibuilder.compile(verbose=False)
