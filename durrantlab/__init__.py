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
        cls.addMGLToolsPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addAutoGrowPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addDeepFragPackage(env, default=bool(cls.getCondaActivationCmd()))

    @classmethod
    def addDeepFragPackage(cls, env, default=False):
        DFRAG_INSTALLED = 'deepFrag_installed'
        dFragCommands = 'git clone {} && cd deepfrag && '.format(cls.getDeepFragGithub())
        dFragCommands += 'conda create -y -n {} -c fastai -c conda-forge --file requirements.txt prody=1.11 rdkit && '\
            .format(cls.getEnvName(DFRAG_DIC))
        dFragCommands += '{} && pip install pyparsing==2.4.7 torch && '.format(cls.getEnvActivationCommand(DFRAG_DIC))
        dFragCommands += "mkdir .store && wget {} -P .store && wget {} -O DFModel.zip && " \
                         "unzip DFModel.zip -d .store/model && ".\
            format(cls.getDeepFragFingerprints(), cls.getDeepFragModel())
        dFragCommands += 'touch ../{}'.format(DFRAG_INSTALLED)
        dFragCommands = [(dFragCommands, DFRAG_INSTALLED)]

        env.addPackage(DFRAG_DIC['name'], version=DFRAG_DIC['version'],
                       tar='void.tgz', commands=dFragCommands, default=default)

    @classmethod
    def addAutoGrowPackage(cls, env, default=False):
        AGROW_INSTALLED = 'autogrow_installed'
        agrowCommands = 'conda create -y -n {} -c rdkit rdkit=2020.09 python=3.7 && '.format(cls.getEnvName(AGROW_DIC))
        agrowCommands += '{} && '.format(cls.getEnvActivationCommand(AGROW_DIC))
        agrowCommands += 'conda install -y numpy=1.21 scipy=1.7 matplotlib=3.5 func_timeout=4.3 && '
        agrowCommands += 'conda install -y -c openbabel openbabel=2.4 && '
        agrowCommands += 'git clone {} && '.format(cls.getAutoGrowGithub())
        agrowCommands += 'touch {}'.format(AGROW_INSTALLED)
        agrowCommands = [(agrowCommands, AGROW_INSTALLED)]

        env.addPackage(AGROW_DIC['name'], version=AGROW_DIC['version'],
                       tar='void.tgz', commands=agrowCommands, default=default)

    @classmethod
    def addMGLToolsPackage(cls, env, default=True):
        # Instantiating install helper
        installer = InstallHelper(MGL_DIC['name'], packageHome=cls.getVar(MGL_DIC['home']),
                                  packageVersion=MGL_DIC['version'])

        mglEnvName = cls.getEnvName(MGL_DIC)

        installer.addCommand(
            f'conda create -y -n {mglEnvName} -c conda-forge -c bioconda mgltools={MGL_DIC["version"]}',
            'MGLTOOLS_ENV_CREATED'
        ).addCommand(
            f'{cls.getEnvActivationCommand(MGL_DIC)} && '
            f'rm -rf {cls.getVar(MGL_DIC["home"])} && '
            f'ln -s $CONDA_PREFIX {cls.getVar(MGL_DIC["home"])}',
            'MGLTOOLS_SYMLINK_CREATED'
        ).addPackage(env, dependencies=['conda'], default=default)

    # ---------------------------------- Utils functions  -----------------------
    @classmethod
    def getPluginHome(cls, path=""):
        import durrantlab
        fnDir = os.path.split(durrantlab.__file__)[0]
        return os.path.join(fnDir, path)

    @classmethod
    def getAutoGrowGithub(cls):
      return 'https://github.com/durrantlab/autogrow4.git'

    @classmethod
    def getDeepFragGithub(cls):
        return 'https://github.com/durrantlab/deepfrag.git'

    @classmethod
    def getDeepFragFingerprints(cls):
        return 'https://durrantlab.pitt.edu/apps/deepfrag/files/fingerprints.h5'

    @classmethod
    def getDeepFragModel(cls):
        return 'https://durrantlab.pitt.edu/apps/deepfrag/files/final_model_v2.zip'

    @classmethod
    def getMGLToolsURL(cls):
        return 'https://ccsb.scripps.edu/download/548/'

