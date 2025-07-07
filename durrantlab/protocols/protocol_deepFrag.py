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
    """DeepFrag lead optimization
    
    User IA Manual: DeepFrag Protocol

The DeepFrag protocol integrates the DeepFrag neural network into the
Scipion-Chem environment, enabling fragment-based ligand optimization driven by
learned structure?activity relationships. DeepFrag is designed to suggest small
chemical modifications to a ligand fragment bound in a protein pocket, using a
grid-based representation of the binding site as input to a trained deep
learning model.

To initiate the protocol, the user must provide a protein?ligand complex in
PDBQT format. This structure should reflect a plausible binding mode, where the
ligand is already positioned within the target receptor. DeepFrag operates by
identifying a modifiable atom within the ligand and proposing new substituents
that may improve binding affinity or introduce structural diversity. Therefore,
one key input is the specification of the ligand atom to be replaced. This atom
must be part of the ligand and should be adjacent to at least one hydrogen to
allow replacement by a valid substituent.

The user defines the grid box parameters that enclose the local environment
around the modification site. This box must be large enough to capture the
chemical context influencing fragment prioritization, including nearby residues
and ligand features. The protocol supports manual definition of the box center
and size, as well as options to derive these values automatically based on the
ligand atom to be replaced.

The model used by DeepFrag is fixed and trained to generate SMILES
representations of suggested fragments. The user can set how many suggestions to
generate per modification site and whether to filter out chemically invalid
fragments based on valence rules or sanitization criteria. Once candidate
fragments are generated, they are attached to the ligand core and saved as new
molecules in standard chemical formats such as SDF or PDBQT.

Output includes the optimized ligand structures with newly attached fragments,
along with a table containing the suggested SMILES strings and associated model
scores. These modified ligands can be directly docked or rescored using other
protocols within Scipion-Chem. In addition, the user can inspect the location
and identity of the modified atoms to assess the structural impact of the
suggestions.

Overall, this protocol enables structure-aware fragment elaboration guided by
deep learning. It supports workflows where an initial ligand is iteratively
refined by generating and evaluating novel chemical variants at specific
positions within the binding site.
    """
    _label = 'DeepFrag lead optimization'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
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
      cId = self._insertFunctionStep('convertStep', prerequisites=[])

      dockSteps = []
      for it, ligandStr in enumerate(self.ligandList.get().split('\n')):
          if ligandStr.strip():
              dockId = self._insertFunctionStep('dFragStep', ligandStr, prerequisites=[cId])
              dockSteps.append(dockId)

      self._insertFunctionStep('createOutputStep', prerequisites=dockSteps)

    def convertStep(self):
        receptorFile = self.getOriginalReceptorFile()
        if receptorFile.endswith('.pdb'):
          shutil.copy(receptorFile, self.getReceptorPDB())
        elif receptorFile.endswith('.pdbqt'):
          self.convertReceptor2PDB(receptorFile)

        oDir = self.getInputLigandsPath()
        if not os.path.exists(oDir):
          os.makedirs(oDir)

        allMols = self.getLigandsFileNames()
        for molFile in allMols:
            inName, inExt = os.path.splitext(os.path.basename(molFile))
            oFile = os.path.join(oDir, inName + '.pdb')

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
            if convMol.endswith('.pdb') and basename in convMol:
              molName, i = convMol.split('__')[0], int(convMol.split('_')[-5])
              oFile = os.path.join(oDir, molName+'.pdb')
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
