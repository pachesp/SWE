# @file
# This file is part of SWE.
#
# @author Kostas Michalopoulos, Andrew Helwer, Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
#
# @section LICENSE
#
# SWE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SWE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SWE.  If not, see <http://www.gnu.org/licenses/>.
#
#
# @section DESCRIPTION
#
# Tool for the usage of CUDA in scons.
#

"""
SCons.Tool.cuda

CUDA Tool for SCons

***
breuera - changes (see original: http://www.scons.org/wiki/CudaTool)
 * commented out *.linkinfo
 * lib64 instead of lib
 * removed unnecessary SDK, "cudaSDKPath + '/lib64', cudaSDKPath + '/common/lib64' + cudaSDKSubLibDir, ", from LIBSPATH
rettenbs
 * removed SDK path completely
***

"""

import os
import sys
import SCons.Tool
import SCons.Scanner.C
import SCons.Defaults

CUDAScanner = SCons.Scanner.C.CScanner()

# this object emitters add '.linkinfo' suffixed files as extra targets
# these files are generated by nvcc. The reason to do this is to remove
# the extra .linkinfo files when calling scons -c
def CUDANVCCStaticObjectEmitter(target, source, env):
        tgt, src = SCons.Defaults.StaticObjectEmitter(target, source, env)
        #for file in src:
                #lifile = os.path.splitext(src[0].rstr())[0] + '.linkinfo'
                #tgt.append(lifile)
        return tgt, src
def CUDANVCCSharedObjectEmitter(target, source, env):
        tgt, src = SCons.Defaults.SharedObjectEmitter(target, source, env)
        #for file in src:
                #lifile = os.path.splitext(src[0].rstr())[0] + '.linkinfo'
                #tgt.append(lifile)
        return tgt, src

def generate(env):
        staticObjBuilder, sharedObjBuilder = SCons.Tool.createObjBuilders(env);
        staticObjBuilder.add_action('.cu', '$STATICNVCCCMD')
        staticObjBuilder.add_emitter('.cu', CUDANVCCStaticObjectEmitter)
        sharedObjBuilder.add_action('.cu', '$SHAREDNVCCCMD')
        sharedObjBuilder.add_emitter('.cu', CUDANVCCSharedObjectEmitter)
        SCons.Tool.SourceFileScanner.add_scanner('.cu', CUDAScanner)

        # default compiler
        env['NVCC'] = 'nvcc'

        # default flags for the NVCC compiler
        env['NVCCFLAGS'] = ''
        env['STATICNVCCFLAGS'] = ''
        env['SHAREDNVCCFLAGS'] = ''
        env['ENABLESHAREDNVCCFLAG'] = '-shared'

        # default NVCC commands
        env['STATICNVCCCMD'] = '$NVCC -o $TARGET -c $NVCCFLAGS $STATICNVCCFLAGS $SOURCES'
        env['SHAREDNVCCCMD'] = '$NVCC -o $TARGET -c $NVCCFLAGS $SHAREDNVCCFLAGS $ENABLESHAREDNVCCFLAG $SOURCES'

        # helpers
        home=os.environ.get('HOME', '')
        programfiles=os.environ.get('PROGRAMFILES', '')
        homedrive=os.environ.get('HOMEDRIVE', '')

        # find CUDA Toolkit path and set CUDA_TOOLKIT_PATH
        try:
                cudaToolkitPath = env['CUDA_TOOLKIT_PATH']
        except:
                paths=[home + '/NVIDIA_CUDA_TOOLKIT',
                       home + '/Apps/NVIDIA_CUDA_TOOLKIT',
                           home + '/Apps/NVIDIA_CUDA_TOOLKIT',
                           home + '/Apps/CudaToolkit',
                           home + '/Apps/CudaTK',
                           '/usr/local/NVIDIA_CUDA_TOOLKIT',
                           '/usr/local/CUDA_TOOLKIT',
                           '/usr/local/cuda_toolkit',
                           '/usr/local/CUDA',
                           '/usr/local/cuda',
                           '/Developer/NVIDIA CUDA TOOLKIT',
                           '/Developer/CUDA TOOLKIT',
                           '/Developer/CUDA',
                           programfiles + 'NVIDIA Corporation/NVIDIA CUDA TOOLKIT',
                           programfiles + 'NVIDIA Corporation/NVIDIA CUDA',
                           programfiles + 'NVIDIA Corporation/CUDA TOOLKIT',
                           programfiles + 'NVIDIA Corporation/CUDA',
                           programfiles + 'NVIDIA/NVIDIA CUDA TOOLKIT',
                           programfiles + 'NVIDIA/NVIDIA CUDA',
                           programfiles + 'NVIDIA/CUDA TOOLKIT',
                           programfiles + 'NVIDIA/CUDA',
                           programfiles + 'CUDA TOOLKIT',
                           programfiles + 'CUDA',
                           homedrive + '/CUDA TOOLKIT',
                           homedrive + '/CUDA']
                pathFound = False
                for path in paths:
                        if os.path.isdir(path):
                                pathFound = True
                                print('scons: CUDA Toolkit found in ' + path)
                                cudaToolkitPath = path
                                break
                if not pathFound:
                        sys.exit("Cannot find the CUDA Toolkit path. Please modify your SConscript or add the path in cudaenv.py")
        env['CUDA_TOOLKIT_PATH'] = cudaToolkitPath

        # add nvcc to PATH
        env.PrependENVPath('PATH', cudaToolkitPath + '/bin')

        # add required libraries
        env.Append(CPPPATH=[cudaToolkitPath + '/include'])
        env.Append(LIBPATH=[cudaToolkitPath + '/lib64'])
        env.Append(RPATH=[cudaToolkitPath + '/lib64'])
        env.Append(LIBS=['cudart'])

def exists(env):
        return env.Detect('nvcc')
