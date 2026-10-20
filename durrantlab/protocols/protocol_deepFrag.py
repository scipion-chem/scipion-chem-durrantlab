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

import os, shutil, json

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, IntParam, TextParam, STEPS_PARALLEL, LabelParam, \
  LEVEL_ADVANCED, StringParam, EnumParam
from pyworkflow.protocol import params
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import runOpenBabel
from pwchem import Plugin as pwchemPlugin
from pwchem.constants import OPENBABEL_DIC, RDKIT_DIC

from durrantlab import Plugin
from durrantlab import DFRAG_DIC, AGROW_DIC
from durrantlab.protocols import ProtChemAutoGrow4

NONE, GYPSUM, OBABEL, RDKIT = 0, 1, 2, 3

class ProtChemDeepFrag(ProtChemAutoGrow4):
    """
    AI Generated:

This protocol runs DeepFrag, a deep convolutional neural network that suggests
fragments to grow a docked ligand from a chosen atom, based on the 3D
environment of the receptor pocket.

The suggested molecules are returned as a SetOfSmallMolecules with their
DeepFrag score, so they can be prepared, docked and compared with the
starting ligand.

Overview
--------
The protocol takes a docked SetOfSmallMolecules (the receptor is taken from
the set) and a list of ligand / connection point pairs:

- Connection point: atom of the ligand where DeepFrag adds a fragment
- Removal point (optional): atom whose fragment is removed first, so that
  DeepFrag proposes a replacement for it

The list is built with the wizards: select the ligand, view it in PyMOL with
its atom names, pick the connection (and removal) atom and add the pair to
the list. Each line of the list is processed independently and in parallel.

Fragment Prediction
-------------------
For each ligand / connection point:

- The receptor and the ligand around the connection point are voxelized
- The network predicts the fingerprint of the ideal fragment for that
  environment, averaging over several grid rotations
- The fragments of the DeepFrag library are ranked by similarity to that
  fingerprint, and the best ones are attached to the ligand

The number of rotations makes the prediction more stable, and the number of
predictions sets how many molecules are returned per pair. DeepFrag runs on
the GPU when enabled, or on the CPU otherwise.

3D Conversion
-------------
DeepFrag returns the full molecules as SMILES. They can be kept as SMILES, or
converted to 3D PDB models with Gypsum-DL (one variant per molecule),
OpenBabel or RDKit.

Workflow
--------
1. Convert the receptor and the listed ligands to PDB (reusing the relabeled
   ligand PDB from the view wizard, so that atom names match)
2. Run DeepFrag for every ligand / connection point
3. Read the predicted molecules and their scores
4. Convert them to 3D with the selected software
5. Register them as SmallMolecules of the output set

Output
------
- outputSmallMolecules:
    One SmallMolecule per suggested molecule, with:

    - _score:
        DeepFrag score of the suggested fragment. Higher values indicate
        a better fit to the pocket environment.
    - molClass:
        "Deepfrag".

Key Features
------------
- Structure-based fragment growing and replacement with deep learning
- Interactive selection of the connection and removal atoms
- Several ligands and growing points in one run, in parallel
- GPU or CPU execution
- SMILES or 3D output (Gypsum-DL, OpenBabel or RDKit)

Use Cases
---------
- Optimizing a docked hit by growing it towards unexplored pocket regions
- Replacing a fragment of a ligand by alternatives better suited to the pocket
- Generating focused analogue libraries for a known binding mode

Notes
-----
- The input molecules must be docked: DeepFrag needs the ligand pose in the
  receptor
- The suggested molecules are not posed in the pocket: their 3D models are
  built from scratch, so they must be docked to compare them with the
  starting ligand
- Atom names must be unique in the ligand: use the view wizard to check them
"""
    _label = 'DeepFrag lead optimization'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
                       label="Use GPU for execution: ",
                       help="This protocol has both CPU and GPU implementation.\
                                                                 Select the one you want to use.")
        form.addHidden(params.GPU_LIST, params.StringParam, default='0', label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used")

        form.addSection(label='Input')
        inGroup = form.addGroup('Input')
        inGroup.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                         label='Input small molecules: ',
                         help="Input small molecules to be docked with AutoGrow")
        inGroup.addParam('inputLigand', StringParam, label='Ligand to optimize: ', default='',
                         help='Specific ligand to optimize using DeepFrag')
        inGroup.addParam('viewLigand', LabelParam, label='View ligand: ',
                         help='Visualize the ligand on PyMol to check the labels')

        inGroup = form.addGroup('Optimization points')
        inGroup.addParam('conPoint', StringParam, label='Connection point: ', default='',
                         help='Specific atom in the ligand where DeepFrag will try to grow new fragments')
        inGroup.addParam('remPoint', StringParam, label='Fragment removal point: ',
                         expertLevel=LEVEL_ADVANCED, default='',
                         help='Specific atom in the ligand whose fragment will be removed from the ligand for fragment '
                              'replacement')

        inGroup.addParam('addLigand', LabelParam, label='Add Ligand / Connection point: ',
                         help='Ligand / Connection point( / Removal point) to the list')
        inGroup.addParam('ligandList', TextParam, width=70,
                         label='List of ligand-connection points:', default='',
                         help='List of ligand-connection points(-remove points) to execute using DeepFrag')

        parGroup = form.addGroup('DeepFrag parameters')
        parGroup.addParam('nGrids', IntParam, label='Grid rotations: ', default=4,
                          help='Number of grid rotations to use. Using more will take longer but produce a more stable '
                               'prediction')
        parGroup.addParam('topK', IntParam, label='Number of predictions: ', default=25,
                          help='Number of predictions to generate per ligand')
        parGroup.addParam('toPdb', EnumParam, label='Convert SMI output to PDB using: ', default=1,
                          choices=['None', 'Gypsum-DL', 'OpenBabel', 'RDKit'],
                          help='Convert SMI output to PDB using the chosen software')

        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
      cId = self._insertFunctionStep(self.convertStep, prerequisites=[])

      dockSteps = []
      for it, ligandStr in enumerate(self.ligandList.get().split('\n')):
          if ligandStr.strip():
              dockId = self._insertFunctionStep(self.dFragStep, ligandStr, prerequisites=[cId])
              dockSteps.append(dockId)

      self._insertFunctionStep(self.createOutputStep, prerequisites=dockSteps)

    def convertStep(self):
        receptorFile = self.getOriginalReceptorFile()
        if receptorFile.endswith('.pdb'):
          shutil.copy(receptorFile, self.getReceptorPDB())
        else:
          self.convertReceptor2PDB(receptorFile)

        oDir = self.getInputLigandsPath()
        if not os.path.exists(oDir):
          os.makedirs(oDir)

        allMols = self.getLigandsFileNames()
        for molFile in allMols:
            inName, inExt = os.path.splitext(os.path.basename(molFile))
            oFile = os.path.join(oDir, inName + '.pdb')
            # Reuse the relabeled PDB if the wizard already produced one
            tmpPdbPath = self.getMolTmpPDBPath(molFile)
            if os.path.exists(tmpPdbPath):
                shutil.copy(tmpPdbPath, oFile)
            else:
                args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(molFile), oFile)
                runOpenBabel(protocol=self, args=args, cwd=oDir)

    def dFragStep(self, ligStr):
        '''Executes AutoGrow to dock and generate ligand variants to look for the best hits'''
        ligDic = self.parseLigandLine(ligStr)
        mol = self.getLigand(ligDic['Ligand'])
        molFile = self.getInputMolPath(mol)
        outName, outDir = '{}_{}'.format(mol.getUniqueName(), ligDic['ConnectionPoint']), self._getExtraPath()
        if 'RemovalPoint' in ligDic:
          outName += '{}'.format(ligDic['RemovalPoint'])

        args = ' --receptor {} --ligand {} --out {}.csv --full --cname {}'. \
          format(self.getReceptorPDB(), molFile, outName, ligDic['ConnectionPoint'])
        if 'RemovalPoint' in ligDic:
          args += ' --rname {}'.format(ligDic['RemovalPoint'])
        args += ' --num_grids {} --top_k {}'.format(self.nGrids.get(), self.topK.get())

        if self.useGpu.get():
            os.environ["CUDA_VISIBLE_DEVICES"] = self.gpuList.get()
        else:
            args += ' --cpu'

        Plugin.runScript(self, 'deepfrag.py', args, env=DFRAG_DIC, cwd=self._getExtraPath(),
                         scriptDir=Plugin.getProgramHome(DFRAG_DIC, path='deepfrag'))

    def createOutputStep(self):
      outDir = os.path.abspath(self._getPath('outputLigands'))
      makePath(outDir)
      outputSet = SetOfSmallMolecules().create(outputPath=outDir)

      for outFile in os.listdir(self._getExtraPath()):
        if outFile.endswith('.csv'):
          molName = outFile.split('.')[0]
          smiDic = self.parseOutMols(self._getExtraPath(outFile))
          print(smiDic)
          molDic = self.convertSMIToPDB(smiDic, outDir, molName)

          for smi in molDic:
            newSmallMol = SmallMolecule(smallMolFilename=molDic[smi][0], molName='guess')
            newSmallMol._score = pwobj.Float(molDic[smi][1])
            newSmallMol.setMolClass('Deepfrag')
            outputSet.append(newSmallMol)

      self._defineOutputs(outputSmallMolecules=outputSet)

    ################################################
    def getOriginalReceptorFile(self):
      return self.inputSmallMolecules.get().getProteinFile()

    def getMolTmpPDBPath(self, molFile):
        '''Project tem path for the relabeled ligand PDB produced by the view wizard'''
        base = os.path.splitext(os.path.basename(molFile))[0]
        return os.path.abspath(self.getProject().getTmpPath('{}.pdb'.format(base)))

    def parseLigandLine(self, ligandLine):
      return json.loads(ligandLine)

    def getLigandsFileNames(self, pose=True):
      ligFns = []
      for mol in self.inputSmallMolecules.get():
        if mol.__str__() in self.ligandList.get():
          if pose and mol.getPoseFile():
            ligFns.append(os.path.abspath(mol.getPoseFile()))
          else:
            ligFns.append(os.path.abspath(mol.getFileName()))
      return ligFns

    def getLigand(self, ligandStr):
      for mol in self.inputSmallMolecules.get():
        if mol.__str__() == ligandStr:
          return mol

    def getInputMolPath(self, mol):
      molFile = self.getInputLigandsPath(os.path.splitext(os.path.basename(mol.getPoseFile()))[0])
      return molFile + '.pdb'

    def parseOutMols(self, csvFile):
      molDic = {}
      with open(csvFile) as f:
        f.readline()
        for line in f:
          molDic[line.split(',')[1]] = float(line.split(',')[-1])
      return molDic

    def convertSMIToPDB(self, smiDic, oDir, basename):
      '''Return {smi: [ConvFile, score]}'''
      molDic = {}
      conv = self.toPdb.get()
      if conv == NONE:
        for i, smi in enumerate(smiDic):
          oFile = os.path.join(oDir, basename + str(i+1) + '.smi')
          with open(oFile, 'w') as f:
            f.write(smi)
          molDic[smi] = [oFile, smiDic[smi]]
      else:
        iDic = self.writeMolsSmiFile(smiDic, basename)

        if conv == GYPSUM:
          args = ' --source {} --output_folder {} --num_processors {} --add_pdb_output --max_variants_per_compound 1'. \
            format(self.getInputSMIFile(), oDir, self.numberOfThreads.get())
          dirPath = 'autogrow4/autogrow/operators/convert_files/gypsum_dl'
          Plugin.runScript(self, 'run_gypsum_dl.py', args, env=AGROW_DIC, cwd=oDir,
                           scriptDir=Plugin.getProgramHome(AGROW_DIC, path=dirPath))

          for convMol in os.listdir(oDir):
              if convMol.endswith('.pdb') and basename and convMol.startswith(basename):
                  suffix = convMol[len(basename):]
                  i = int(suffix.lstrip('_').split('_')[0])
                  oFile = os.path.join(oDir, '{}_{}.pdb'.format(basename, i))
                  os.rename(os.path.join(oDir, convMol), oFile)
                  molDic[iDic[i]] = [oFile, smiDic[iDic[i]]]

        else:
          convScript = 'obabel_IO.py' if conv == OBABEL else 'rdkit_IO.py'
          convEnv = OPENBABEL_DIC if conv == OBABEL else RDKIT_DIC

          args = ' -i "{}" -of pdb --outputDir "{}"'.format(self.getInputSMIFile(), oDir)
          args += ' --make3D -nt {}'.format(self.numberOfThreads.get())
          pwchemPlugin.runScript(self, convScript, args, env=convEnv, cwd=oDir)

          for oFile in os.listdir(oDir):
            if oFile.endswith('.pdb'):
              i = int(oFile.split('_')[-1].split('.')[0])
              molDic[iDic[i]] = [os.path.join(oDir, oFile), smiDic[iDic[i]]]

      return molDic

    def writeMolsSmiFile(self, smiDic, basename):
      with open(self.getInputSMIFile(), 'w') as f:
        iDic = {}
        for i, smi in enumerate(smiDic):
          f.write('{}\t{}_{}\n'.format(smi, basename, i + 1))
          iDic[i + 1] = smi
      return iDic

    def _validate(self):
      vals = []
      if not self.inputSmallMolecules.get().isDocked():
        vals.append('The Input small molecules must be docked to a protein for using DeepFrag')

      if not self.ligandList.get().strip():
        vals.append('You must add some ligands-connection points to the "Ligand List" input. To do so, use '
                    'the wizards to select the ligand and the connection (and removal) points and add it with the '
                    '"Add ligand" wizard.')
      return vals

    def _warnings(self):
      warns = []
      return warns
