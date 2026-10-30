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

import os, re, shutil

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, IntParam, FloatParam, STEPS_PARALLEL, BooleanParam, \
  LEVEL_ADVANCED, StringParam, EnumParam
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import runOpenBabel, calculate_centerMass

from durrantlab import Plugin
from durrantlab import AGROW_DIC

argsSkipDic = {'skip_optimize_geometry': 'optGeo', 'skip_alternate_ring_conformations': 'ringConf',
               'skip_adding_hydrogen': 'addHyd', 'skip_making_tautomers': 'makeTau',
               'skip_enumerate_chiral_mol': 'enumChir', 'skip_enumerate_double_bonds': 'enumBonds'}
argsBoolDic = {'let_tautomers_change_chirality': 'chirChang', 'use_durrant_lab_filters': 'durrFilter'}

class ProtChemGypsumDL(EMProtocol):
    """
    AI Generated:

This protocol prepares small molecules with Gypsum-DL: it enumerates their
probable ionization states, tautomers, stereoisomers and ring conformations,
and generates 3D models for each variant.

The variants are returned as a new SetOfSmallMolecules ready to be docked or
used by other virtual screening protocols.

Overview
--------
The protocol takes a SetOfSmallMolecules in any format OpenBabel can read.
All the molecules are converted with OpenBabel into a single SMILES file,
which is the input of Gypsum-DL. Gypsum-DL is run from the copy bundled with
AutoGrow4, inside the AutoGrow4 conda environment.

Since the input is converted to SMILES, its original 3D coordinates are not
used: Gypsum-DL builds new 3D models from scratch.

Variant Enumeration
-------------------
For each input molecule, Gypsum-DL can:

- Protonate it for the selected pH range (pH 6.4 - 8.4 by default, in steps
  defined by the pKa precision)
- Enumerate its tautomers
- Enumerate unspecified chiral centers and double bonds
- Generate alternative non-aromatic ring conformations
- Optimize the geometry of each 3D model

Each of these operations can be switched off. Optionally, tautomers that
change the number of chiral centers can be allowed, and the Durrant lab
substructure filters can remove variants judged improbable.

The number of variants kept per molecule (5 by default) and the
exhaustiveness of the conformer search (3 by default) control the size of
the output and the run time.

Workflow
--------
1. Convert all the input molecules to a single SMILES file (OpenBabel)
2. Run Gypsum-DL on it, using the protocol threads as processors
3. Rename each output PDB as <molName>-<n>.pdb
4. Register each variant as a SmallMolecule of the output set

Output
------
- outputSmallMolecules:
    One SmallMolecule per generated variant, in PDB format, with:

    - molName:
        Name of the input molecule the variant comes from.
    - confId:
        Index of the variant for that molecule (1, 2, ...).
    - type:
        "GypsumDL".

Key Features
------------
- Protonation, tautomer, stereoisomer and ring conformer enumeration in one
  step
- Configurable pH range and pKa precision
- Optional Durrant lab filters for improbable variants
- Parallel execution over the protocol threads

Use Cases
---------
- Preparing ligand libraries before docking
- Exploring the protonation and tautomeric states of a hit
- Generating 3D models for molecules only available as SMILES

Notes
-----
- The molecule name is taken from the text before the first "_" of the
  Gypsum-DL output file, so input names containing "_" are merged together
- Gypsum-DL runs in the AutoGrow4 environment, so AutoGrow4 must be
  installed
"""
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

      self.writeInputSMIFile(self.getLigandsFileNames(), self.getInputSMIFile())

    def writeInputSMIFile(self, molFiles, smiFile, batchSize=300):
      '''Writes the molecules to a SMILES file with the largest fragment of each one (salts and solvents removed)
      and a safe name. OpenBabel names can be file paths (PDB inputs) or SDF fields spanning several lines, which
      break the Gypsum-DL / AutoGrow4 file handling'''
      usedNames = set()
      with open(smiFile, 'w') as f:
        for i in range(0, len(molFiles), batchSize):
          batch = molFiles[i:i + batchSize]
          tmpFile = '{}_batch_{}.smi'.format(smiFile, i)
          runOpenBabel(self, '{} -osmi -r -xn -O {}'.format(' '.join(batch), tmpFile), popen=True)
          smis = []
          if os.path.exists(tmpFile):
            with open(tmpFile) as tf:
              smis = [line.split()[0] for line in tf if line.strip()]
            os.remove(tmpFile)

          # Each SMILES is named after its file, unless some file held several molecules
          if len(smis) == len(batch):
            names = [re.sub(r'[^A-Za-z0-9-]', '-', os.path.splitext(os.path.basename(fn))[0]) for fn in batch]
          else:
            names = ['mol{}'.format(i + j + 1) for j in range(len(smis))]

          for smi, name in zip(smis, names):
            uniqueName, k = name, 1
            while uniqueName in usedNames:
              k += 1
              uniqueName = '{}-{}'.format(name, k)
            usedNames.add(uniqueName)
            f.write('{}\t{}\n'.format(smi, uniqueName))

    def getBoolArgs(self):
      args = ''
      for key, pKey in argsSkipDic.items():
        if not getattr(self, pKey).get():
          args += ' --{}'.format(key)

      for key, pKey in argsBoolDic.items():
        if getattr(self, pKey).get():
          args += ' --{}'.format(key)
      return args
