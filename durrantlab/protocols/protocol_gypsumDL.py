# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import os, shutil

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, IntParam, FloatParam, STEPS_PARALLEL, BooleanParam, \
  LEVEL_ADVANCED, StringParam, EnumParam
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import runOpenBabel, calculate_centerMass, getBaseFileName

from durrantlab import Plugin
from durrantlab import AGROW_DIC

argsSkipDic = {'skip_optimize_geometry': 'optGeo', 'skip_alternate_ring_conformations': 'ringConf',
               'skip_adding_hydrogen': 'addHyd', 'skip_making_tautomers': 'makeTau',
               'skip_enumerate_chiral_mol': 'enumChir', 'skip_enumerate_double_bonds': 'enumBonds'}
argsBoolDic = {'let_tautomers_change_chirality': 'chirChang', 'use_durrant_lab_filters': 'durrFilter'}

class ProtChemGypsumDL(EMProtocol):
    """Protocol to execute Gypsum-DL for ligand preparation"""
    _label = 'Gypsum ligand preparation'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineGypsumDLParams(self, form):
      lpGroup = form.addGroup('Ligand preparation')
      lpGroup.addParam('maxVars', IntParam, label='Maximum number of variants per compound: ', default=5,
                       help='Maxmimum number of conformers made per ligand by Gypsum-DL')
      lpGroup.addParam('confExhaust', IntParam, label='Conformer search exhaustiveness: ', default=3,
                       help='How widely Gypsum-DL will search for low-energy conformers. Larger values increase run '
                            'times but can produce better results. Empty for default')
      line = lpGroup.addLine('PH range: ', expertLevel=LEVEL_ADVANCED,
                             help='PH range for Gypsum-DL to consider in the conformer generation. '
                                  'PKa precision: Size of pH substructure ranges.')
      line.addParam('minPH', FloatParam, label='Minimum: ', default=6.4)
      line.addParam('maxPH', FloatParam, label='Maximum: ', default=8.4)
      line.addParam('pkaPrecision', FloatParam, label='PKa presicion: ', default=1.0)
      return lpGroup

    def _defineParams(self, form):
        form.addSection(label='Input')
        inputGroup = form.addGroup('Input small molecules')
        inputGroup.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Input small molecules: ',
                       help="Input small molecules to be docked with AutoGrow")
        prepGroup = self._defineGypsumDLParams(form)
        prepGroup.addParam('optGeo', BooleanParam, label='Optimize geometry: ', default=True,
                           expertLevel=LEVEL_ADVANCED, help='Optimize the small molecule geometry')
        prepGroup.addParam('ringConf', BooleanParam, label='Generate non-aromatic ring conformations: ', default=True,
                           expertLevel=LEVEL_ADVANCED, help='Generate non-aromatic ring conformations')
        prepGroup.addParam('addHyd', BooleanParam, label='Ionizate molecules: ', default=True,
                           expertLevel=LEVEL_ADVANCED, help='Add hydrogens to ionizate the molecules')
        prepGroup.addParam('makeTau', BooleanParam, label='Generate tautomers: ', default=True,
                           expertLevel=LEVEL_ADVANCED, help='Generate tautormers of the small molecules')
        prepGroup.addParam('enumChir', BooleanParam, label='Enumerate chiral centers: ', default=True,
                           expertLevel=LEVEL_ADVANCED, help='Enumeration of unspecified chiral centers')
        prepGroup.addParam('enumBonds', BooleanParam, label='Enumerate double bonds: ', default=True,
                           expertLevel=LEVEL_ADVANCED, help='Enumeration of double bonds')
        prepGroup.addParam('chirChang', BooleanParam, label='Allow change in chiral centers: ',
                           expertLevel=LEVEL_ADVANCED, default=False,
                           help='Allow tautomers that change the total number of chiral centers.')
        prepGroup.addParam('durrFilter', BooleanParam, label='Use Durrant Lab filters: ',
                           expertLevel=LEVEL_ADVANCED, default=False,
                           help='Use substructure filters to remove molecular variants that, though technically '
                                'possible, were judged improbable by members of the Durrant lab.')

        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
      cId = self._insertFunctionStep('convertStep', prerequisites=[])
      prepId = self._insertFunctionStep('preparationStep', prerequisites=[cId])
      self._insertFunctionStep('createOutputStep', prerequisites=[prepId])

    def convertStep(self):
        self.convertInputSMIFile()

    def preparationStep(self):
        '''Executes AutoGrow to dock and generate ligand variants to look for the best hits'''
        outDir = os.path.abspath(self._getPath('outputLigands'))
        if not os.path.exists(outDir):
          os.makedirs(outDir)

        args = ' --source {} --output_folder {} --num_processors {} --add_pdb_output'.\
          format(self.getInputSMIFile(), outDir, self.numberOfThreads.get())
        args += ' --max_variants_per_compound {} --thoroughness {}'.format(self.maxVars.get(), self.confExhaust.get())
        args += ' --min_ph {} --max_ph {} --pka_precision {}'.\
          format(self.minPH.get(), self.maxPH.get(), self.pkaPrecision.get())
        args += self.getBoolArgs()

        dirPath = 'autogrow4/autogrow/operators/convert_files/gypsum_dl'
        Plugin.runScript(self, 'run_gypsum_dl.py', args, env=AGROW_DIC, cwd=outDir,
                         scriptDir=Plugin.getProgramHome(AGROW_DIC, path=dirPath))

    def createOutputStep(self):
      outDir = self._getPath('outputLigands')
      outputSet = SetOfSmallMolecules().create(outputPath=outDir)

      countDic = {}
      for outMol in os.listdir(outDir):
        if outMol.endswith('.pdb'):
          molName = outMol.split('_')[0]
          if molName not in countDic :
            countDic[molName] = 1
          else:
            countDic[molName] += 1

          molFile = os.path.join(outDir, '{}-{}.pdb'.format(molName, countDic[molName]))
          os.rename(os.path.join(outDir, outMol), molFile)

          newSmallMol = SmallMolecule(smallMolFilename=molFile, type='GypsumDL')
          newSmallMol.setMolName(molName)
          newSmallMol.setConfId(countDic[molName])
          outputSet.append(newSmallMol)

      self._defineOutputs(outputSmallMolecules=outputSet)
      self._defineSourceRelation(self.inputSmallMolecules, outputSet)

    ################################################
    # todo: tests
    def getInputLigandsPath(self, path=''):
        return os.path.abspath(self._getExtraPath('inputLigands', path))

    def getInputSMIFile(self):
        return os.path.join(self.getInputLigandsPath(), 'inputSourceCompounds.smi')

    def getLigandsFileNames(self, pose=False):
      ligFns = []
      for mol in self.inputSmallMolecules.get():
        if pose and mol.getPoseFile():
          ligFns.append(os.path.abspath(mol.getPoseFile()))
        else:
          ligFns.append(os.path.abspath(mol.getFileName()))
      return ligFns

    def convertInputSMIFile(self):
      '''Convert source compounds to SMI and store them into a file'''
      oDir = self.getInputLigandsPath()
      if not os.path.exists(oDir):
        os.mkdir(oDir)

      allMols = self.getLigandsFileNames()
      runOpenBabel(self, '{} -O {}'.format(' '.join(allMols), self.getInputSMIFile()),
                   popen=True)

    def getBoolArgs(self):
      args = ''
      for key, pKey in argsSkipDic.items():
        if not getattr(self, pKey).get():
          args += ' --{}'.format(key)

      for key, pKey in argsBoolDic.items():
        if getattr(self, pKey).get():
          args += ' --{}'.format(key)
      return args
