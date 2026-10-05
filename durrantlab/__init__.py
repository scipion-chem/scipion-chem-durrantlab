# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os, subprocess

import pwchem
from scipion.install.funcs import InstallHelper

from .constants import *

_version_ = '0.1'
_logo = "durrant_logo.png"
_references = ['']

class Plugin(pwchem.Plugin):
    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(AGROW_DIC['home'], '{}-{}'.format(AGROW_DIC['name'], AGROW_DIC['version']))
        cls._defineEmVar(DFRAG_DIC['home'], '{}-{}'.format(DFRAG_DIC['name'], DFRAG_DIC['version']))
        cls._defineVar("AUTOGROW_ENV_ACTIVATION", cls.getEnvActivationCommand(AGROW_DIC))
        cls._defineVar("DEEPFRAG_ENV_ACTIVATION", cls.getEnvActivationCommand(DFRAG_DIC))

    @classmethod
    def defineBinaries(cls, env):
        # pwchem's MGLTools: skipped if already declared or already installed in software/em
        pwchem.Plugin.addMGLToolsPackage(env)
        cls.addAutoGrowPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addDeepFragPackage(env, default=bool(cls.getCondaActivationCmd()))

    @classmethod
    def addDeepFragPackage(cls, env, default=False):
        installer = InstallHelper(DFRAG_DIC['name'], packageHome=cls.getVar(DFRAG_DIC['home']),
                                  packageVersion=DFRAG_DIC['version'])

        envName = cls.getEnvName(DFRAG_DIC)
        envActivation = cls.getEnvActivationCommand(DFRAG_DIC)

        # Only the inference deps: deepfrag's requirements.txt is a 2020 training freeze that takes hours to solve
        installer.addCommand(
            f'git clone {cls.getDeepFragGithub()} && '
            f'conda create -y -n {envName} --override-channels -c conda-forge '
            f'python=3.11 "numpy<2" h5py numba=0.61 rdkit=2023.09 prody tqdm requests',
            'DEEPFRAG_ENV_CREATED'
        ).addCommand(
            # numba's built-in CUDA fails on CUDA 13 drivers: use numba-cuda with its ctypes driver binding
            f'{envActivation} && pip install torch==2.6.0 "numba-cuda[cu12]==0.15.2" && '
            f'conda env config vars set NUMBA_CUDA_USE_NVIDIA_BINDING=0',
            'DEEPFRAG_DEPS_INSTALLED'
        ).addCommand(
            f'cd deepfrag && '
            f'mkdir -p .store && '
            f'wget {cls.getDeepFragFingerprints()} -P .store && '
            f'wget {cls.getDeepFragModel()} -O DFModel.zip && '
            f'unzip DFModel.zip -d .store/model && '
            f'rm DFModel.zip',
            'DEEPFRAG_MODELS_DOWNLOADED'
        ).addPackage(env, dependencies=['conda'], default=default)

    @classmethod
    def addAutoGrowPackage(cls, env, default=False):
        installer = InstallHelper(AGROW_DIC['name'], packageHome=cls.getVar(AGROW_DIC['home']),
                                  packageVersion=AGROW_DIC['version'])

        envName = cls.getEnvName(AGROW_DIC)
        envActivation = cls.getEnvActivationCommand(AGROW_DIC)

        installer.addCommand(
            f'conda create -y -n {envName} -c rdkit rdkit=2020.09 python=3.7',
            'AUTOGROW_ENV_CREATED'
        ).addCommand(
            f'{envActivation} && '
            f'conda install -y numpy=1.21 scipy=1.7 matplotlib=3.5 func_timeout=4.3 && '
            f'conda install -y -c openbabel openbabel=2.4',
            'AUTOGROW_DEPS_INSTALLED'
        ).addCommand(
            f'wget {cls.getAutoGrowUrl()} -O autogrow4.zip && '
            f'unzip autogrow4.zip && '
            f'mv autogrow4-4.0.3 autogrow4 && '
            f'rm autogrow4.zip',
            'AUTOGROW_DOWNLOADED'
        ).addPackage(env, dependencies=['conda'], default=default)

    # ---------------------------------- Utils functions  -----------------------
    @classmethod
    def getPluginHome(cls, path=""):
        import durrantlab
        fnDir = os.path.split(durrantlab.__file__)[0]
        return os.path.join(fnDir, path)

    @classmethod
    def getAutoGrowUrl(cls):
        return 'https://github.com/durrantlab/autogrow4/archive/refs/tags/v4.0.3.zip'

    @classmethod
    def getDeepFragGithub(cls):
        return 'https://github.com/durrantlab/deepfrag.git'

    @classmethod
    def getDeepFragFingerprints(cls):
        return 'https://durrantlab.pitt.edu/apps/deepfrag/files/fingerprints.h5'

    @classmethod
    def getDeepFragModel(cls):
        return 'https://durrantlab.pitt.edu/apps/deepfrag/files/final_model_v2.zip'