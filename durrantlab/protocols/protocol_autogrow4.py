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
from durrantlab.protocols.protocol_gypsumDL import ProtChemGypsumDL

argsDic = {'max_variants_per_compound': 'maxVars', # 'gypsum_thoroughness': 'confExhaust',
           'min_ph': 'minPH', 'max_ph': 'maxPH', 'pka_precision': 'pkaPrecision',
           'min_atom_match_MCS': 'minMatch',
           'docking_exhaustiveness': 'dockExhaust', 'docking_num_modes': 'nDocks',
           'num_generations': 'nGens', 'tourn_size': 'tournSize', 'number_of_crossovers': 'nCrossovers',
           'number_of_mutants': 'nMutants', 'number_elitism_advance_from_previous_gen': 'nElitism',
           'top_mols_to_seed_next_generation': 'topMolsSeed'}

argsBoolDic = {'protanate_step': 'protonate',
               'MozziconacciFilter': 'fMozziconacci',
               'VandeWaterbeemdFilter': 'fVandeWaterbeemd', 'PAINSFilter': 'fPAINS', 'NIHFilter': 'fNIH',
               'BRENKFilter': 'fBRENK'}

argsBoolEnumDic = {'filter_source_compounds': 'fSource'}

argsEnumDic = {'rxn_library': 'reactionLibrary', 'dock_choice': 'dockChoice', 'scoring_choice': 'scoreChoice',
               'selector_choice': 'selectorChoice'}
zincDic = {'ZINC MW <= 100 Da': 'Fragment_MW_up_to_100.smi',
               'ZINC 100 Da < MW <= 150 Da': 'Fragment_MW_100_to_150.smi',
               'ZINC 150 Da < MW <= 200 Da': 'Fragment_MW_150_to_200.smi',
               'ZINC 200 Da < MW <= 250 Da': 'Fragment_MW_200_to_250.smi',
               'Naphthalene': 'naphthalene_smiles.smi', 'PARPi': 'PARPi.smi',
               'PARPi_BRICS_frags': 'PARPi_BRICS_frags.smi'}

class ProtChemAutoGrow4(ProtChemGypsumDL):
    """AutoGrow for docking and lead optimization"""
    _label = 'AutoGrow4 docking and lead optimization'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        inputGroup = form.addGroup('Receptor specification')
        inputGroup.addParam('fromReceptor', EnumParam, label='Dock on : ', default=1,
                       choices=['Whole protein', 'SetOfStructROIs'], display=EnumParam.DISPLAY_HLIST,
                       help='Whether to dock on a whole protein surface or on specific regions')

        # Docking on whole protein
        inputGroup.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct",
                       label='Input atomic structure: ', condition='fromReceptor == 0',
                       help="The atom structure to use as receptor in the docking")
        inputGroup.addParam('radius', FloatParam, label='Grid radius for whole protein: ',
                       condition='fromReceptor == 0', allowsNull=False,
                       help='Radius of the Autodock grid for the whole protein')

        # Docking on pockets
        inputGroup.addParam('inputStructROIs', PointerParam, pointerClass="SetOfStructROIs",
                       label='Input pockets: ', condition='not fromReceptor == 0',
                       help="The protein structural ROIs to dock in")
        inputGroup.addParam('pocketRadiusN', FloatParam, label='Grid radius vs StructROI radius: ',
                       condition='not fromReceptor == 0', default=1.1, allowsNull=False,
                       help='The radius * n of each StructROI will be used as grid radius')

        inputGroup.addParam('ligandsSource', EnumParam, label='Ligands origin : ', default=1,
                            choices=['AutoGrow libraries', 'User small molecules'], display=EnumParam.DISPLAY_HLIST,
                            help='Source of the input molecules that will be variated and docked')
        inputGroup.addParam('inputAGrowMols', EnumParam, label='Input AutoGrow small molecules: ',
                           condition='ligandsSource==0',
                           choices=['ZINC MW <= 100 Da', 'ZINC 100 Da < MW <= 150 Da', 'ZINC 150 Da < MW <= 200 Da',
                                    'ZINC 200 Da < MW <= 250 Da', 'Naphthalene', 'PARPi', 'PARPi_BRICS_frags'],
                           help="Source compounds from AutoGrow examples. For more info about the groups, check "
                                "Example_source_compound_notes.txt file in autogrow4/source_compounds folder")
        inputGroup.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Input small molecules: ', condition='ligandsSource==1',
                       help="Input small molecules to be docked with AutoGrow")

        lpGroup = form.addGroup('Ligand preparation')
        lpGroup.addParam('maxVars', IntParam, label='Maximum number of variants per compound: ', default=3,
                         help='Maxmimum number of conformers made per ligand by Gypsum-DL')
        # lpGroup.addParam('confExhaust', StringParam, label='Conformer search exhaustiveness: ', default='',
        #                  help='How widely Gypsum-DL will search for low-energy conformers. Larger values increase run '
        #                       'times but can produce better results. Empty for default')
        line = lpGroup.addLine('PH range: ', expertLevel=LEVEL_ADVANCED,
                              help='PH range for Gypsum-DL to consider in the conformer generation. '
                                   'PKa precision: Size of pH substructure ranges.')
        line.addParam('minPH', FloatParam, label='Minimum: ', default=6.4)
        line.addParam('maxPH', FloatParam, label='Maximum: ', default=8.4)
        line.addParam('pkaPrecision', FloatParam, label='PKa presicion: ', default=1.0)

        form.addSection('Reactions')
        smGroup = form.addGroup('SmilesMerge')
        smGroup.addParam('minMatch', IntParam, label='Minimum atoms for MCS match: ', default=4,
                         help='Determines the minimum number of atoms in common for a substructurematch. The higher '
                              'the more restrictive, but the more likely for two ligands not to match')
        smGroup.addParam('protonate', BooleanParam, label='Protonate step: ', default=False,
                        help='Indicates if Smilesmerge uses protanated mols (if true) or deprot (if False) '
                             'SmilesMerge is 10x faster when deprotanated')
        smGroup = form.addGroup('Mutations')
        smGroup.addParam('reactionLibrary', EnumParam, label='Reactions library : ', default=2,
                         choices=["click_chem_rxns", "robust_rxns", "all_rxns"],
                         help='This set of reactions to be used in Mutation.')

        form.addSection('Docking and Genetic algorithm')
        dGroup = form.addGroup('Docking')
        dGroup.addParam('dockChoice', EnumParam, label='Docking software: ', default=1,
                        choices=["VinaDocking", "QuickVina2Docking"], display=EnumParam.DISPLAY_HLIST,
                        help='Which molecular docking software to use')
        dGroup.addParam('dockExhaust', StringParam, label='Docking exhaustiveness: ',
                        default='', expertLevel=LEVEL_ADVANCED,
                        help='Exhaustiveness of the global search (roughly proportional to time). If empty, default for'
                             ' each software will be used')
        dGroup.addParam('nDocks', StringParam, label='Number of docking modes: ', default='',
                        help='Maximum number of binding modes to generate in docking. If empty, default for'
                             ' each software will be used')
        dGroup.addParam('scoreChoice', EnumParam, label='Scoring choice: ', default=0,
                        choices=["VINA", "NN1", "NN2"], expertLevel=LEVEL_ADVANCED,
                        help='The scoring_choice to use to assess the ligands docking fitness. Default is using '
                             'Vina/QuickVina2 ligand affinity while NN1/NN2 use a Neural Network to assess the docking '
                             'pose')

        gaGroup = form.addGroup('Genetic algorithm')
        gaGroup.addParam('nGens', IntParam, label='Number of generations: ', default=10,
                          help='Number of generations of different ligand variants for AutoGrow to generate')
        gaGroup.addParam('genMin', IntParam, label='Minimum generation for output: ', default=1,
                         expertLevel=LEVEL_ADVANCED,
                         help='Minimum number of generations to be performed before considering an output docking pose')
        gaGroup.addParam('selectorChoice', EnumParam, label='Selector choice : ', default=0,
                         choices=["Roulette_Selector", "Rank_Selector", "Tournament_Selector"],
                         help='This determines whether the fitness criteria are chosen by a Weighted Roulette, Ranked, '
                              'or Tournament style Selector. The Rank option is a non-redundant selector. Roulette and '
                              'Tournament chose without replacement and are stoichastic options. Warning do not use '
                              'Rank_Selector for small runs as there is potential that the number of desired ligands '
                              'exceed the number of ligands to chose from."')
        gaGroup.addParam('tournSize', FloatParam, label='Tournament size: ', default=0.1,
                         condition='selectorChoice == 2', expertLevel=LEVEL_ADVANCED,
                         help='If using the Tournament_Selector this determines the size of each tournament. The '
                              'number of ligands used for each tournament will the tourn_size * the number of '
                              'considered ligands.')
        gaGroup.addParam('nCrossovers', IntParam, label='Number of crossovers: ',
                         default=10, expertLevel=LEVEL_ADVANCED,
                         help='The number of ligands which will be created via crossover in each generation')
        gaGroup.addParam('nMutants', IntParam, label='Number of mutants: ',
                         default=10, expertLevel=LEVEL_ADVANCED,
                         help='The number of ligands which will be created via mutation in each generation')
        gaGroup.addParam('nElitism', IntParam, label='Number of elitists: ',
                         default=10, expertLevel=LEVEL_ADVANCED,
                         help='The number of ligands chosen for elitism. These will advance from the previous '
                              'generation directly into the next generation. This is purely advancing based on '
                              'Docking/Rescore fitness. This does not select for diversity.')

        gaGroup.addParam('topMolsSeed', IntParam, label='Top mols to seed: ',
                         default=10, expertLevel=LEVEL_ADVANCED,
                         help='Number of mols that seed next generation.')

        form.addSection('Filtering')
        form.addParam('fSource', BooleanParam, label='Filter source ligands: ', default=True,
                      help='Whether to use the filters defined in this form on the source ligands or just in those '
                      'created by the genetic algorithm in subsequent generations')

        fGroup = form.addGroup('Molecular filters')
        fGroup.addParam('fLipinski', BooleanParam, label='Lipinski filter: ', default=False,
                        help='Lipinski filters for orally available drugs following Lipinski rule of fives. Filters by '
                             'molecular weight, logP and number of hydrogen bond donors and acceptors. Strict '
                             'implementation means a ligand must pass all requirements. Lenient implementation means a '
                             'ligand may fail all but one requirement and still passes.')
        fGroup.addParam('LipinskiType', EnumParam, label='Lipinski filter type: ', default=0,
                            choices=['Strict', 'Lenient'], display=EnumParam.DISPLAY_HLIST, condition='fLipinski',
                            help='Type of Lipinski filter to apply')
        fGroup.addParam('fGhose', BooleanParam, label='Ghose filter: ', default=False,
                        help='Ghose filters for drug-likeliness; filters by molecular weight, logP and number of atoms')
        fGroup.addParam('GhoseType', EnumParam, label='Ghose filter type: ', default=0,
                        choices=['Regular', 'Modified'], display=EnumParam.DISPLAY_HLIST, condition='fGhose',
                        help='Type of Ghose filter to apply')
        fGroup.addParam('fMozziconacci', BooleanParam, label='Mozziconacci filter: ', default=False,
                        help='Mozziconacci filters for drug-likeliness; filters by the number of rotatable bonds, '
                             'rings, oxygens, and halogens.')
        fGroup.addParam('fVandeWaterbeemd', BooleanParam, label='VandeWaterbeemd filter: ',
                        default=False, expertLevel=LEVEL_ADVANCED,
                        help='VandeWaterbeemd filters for drug likely to be blood brain barrier permeable. Filters by '
                             'the number of molecular weight and Polar Surface Area (PSA).')
        fGroup.addParam('fPAINS', BooleanParam, label='PAINS filter: ', default=False, expertLevel=LEVEL_ADVANCED,
                        help='PAINS filters against Pan Assay Interference Compounds using substructure a search.')
        fGroup.addParam('fNIH', BooleanParam, label='NIH filter: ', default=False, expertLevel=LEVEL_ADVANCED,
                        help='NIH filters against molecules with undersirable functional groups using substructure '
                             'a search.')
        fGroup.addParam('fBRENK', BooleanParam, label='BRENK filter: ', default=False, expertLevel=LEVEL_ADVANCED,
                        help='BRENK filter for lead-likeliness, by matching common false positive molecules to the '
                             'current mol.')

        fGroup.addParam('maxEnergy', FloatParam, label='Maximum energy for docking pose: ',
                        default=0,
                        help='Maximum energy for a docking pose to be considered as output')

        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
      cId = self._insertFunctionStep('convertStep', prerequisites=[])

      dockSteps = []
      if self.fromReceptor.get() == 0:
        dockId = self._insertFunctionStep('dockStep', prerequisites=[cId])
        dockSteps.append(dockId)
      else:
        for it, pocket in enumerate(self.inputStructROIs.get()):
          dockId = self._insertFunctionStep('dockStep', pocket.clone(), it=it, prerequisites=[cId])
          dockSteps.append(dockId)

      self._insertFunctionStep('createOutputStep', prerequisites=dockSteps)

    def convertStep(self):
        self.convertInputSMIFile()

        receptorFile = self.getOriginalReceptorFile()
        if receptorFile.endswith('.pdb'):
          shutil.copy(receptorFile, self.getReceptorPDB())
        elif receptorFile.endswith('.pdbqt'):
          self.convertReceptor2PDB(receptorFile)

    def dockStep(self, pocket=None):
        '''Executes AutoGrow to dock and generate ligand variants to look for the best hits'''
        outDir = self.getOutputPocketDir(pocket)
        if not os.path.exists(outDir):
          makePath(outDir)

        if self.fromReceptor.get() == 0:
          radius = self.radius.get()
          structure, x_center, y_center, z_center = calculate_centerMass(self.getReceptorPDB())
          nCPUs = self.numberOfThreads.get()
        else:
          radius = (pocket.getDiameter() / 2) * self.pocketRadiusN.get()
          x_center, y_center, z_center = pocket.calculateMassCenter()
          nCPUs = self.numberOfThreads.get() // len(self.inputStructROIs.get())

        args = ' --filename_of_receptor {} --source_compound_file {} --root_output_folder {} -p {}'. \
          format(self.getReceptorPDB(), self.getInputSMIFile(), outDir, nCPUs)
        args += ' --center_x {} --center_y {} --center_z {} --size_x {} --size_y {} --size_z {} '.\
          format(x_center, y_center, z_center, radius, radius, radius)
        args += self.getConversionArgs()
        args += self.getArgs()

        Plugin.runScript(self, 'RunAutogrow.py', args, env='AutoGrow4',
                                  cwd=self._getExtraPath(),
                                  scriptDir=Plugin.getProgramHome(AGROW_DIC, path='autogrow4'))

    def createOutputStep(self):
      outDir = self._getPath('outputLigands')
      makePath(outDir)
      outputSet = SetOfSmallMolecules().create(outputPath=outDir)

      for pocketDir in self.getPocketDirs():
        gridId = self.getGridId(pocketDir)
        for genDir in self.getGenDirs(pocketDir, self.genMin.get()):
          dockDic = self.getDockedDic(genDir)
          self.compressResults(genDir)

          for model in dockDic:
            energy = self.parseModelEnergy(dockDic[model])
            if energy <= self.maxEnergy.get():
              posFile = os.path.join(outDir, model + '.pdbqt')
              with open(posFile, 'w') as f:
                f.write(dockDic[model])

              newSmallMol = SmallMolecule(smallMolFilename=posFile, molName='guess')
              newSmallMol._energy = pwobj.Float(energy)

              posId = model.split('_')[-1]
              newSmallMol.poseFile.set(posFile)
              newSmallMol.setPoseId(posId)

              newSmallMol.gridId.set(gridId)
              newSmallMol.setMolClass('AutoGrow4')
              newSmallMol.setDockId(self.getObjId())

              outputSet.append(newSmallMol)

      outputSet.proteinFile.set(self.getOriginalReceptorFile())
      outputSet.setDocked(True)
      self._defineOutputs(outputSmallMolecules=outputSet)
      if self.ligandsSource.get() == 1:
        self._defineSourceRelation(self.inputSmallMolecules, outputSet)


    ################################################

    def getArgs(self):
        args = ' --reduce_files_sizes False'
        for key, protKey in argsDic.items():
          argVal = getattr(self, protKey).get()
          if argVal:
            argVal = int(float(argVal)) if key == 'gypsum_thoroughness' else argVal
            args += ' --{} {}'.format(key, argVal)

        for key, protKey in argsBoolDic.items():
          if getattr(self, protKey).get():
            args += ' --{}'.format(key)

        for key, protKey in argsBoolEnumDic.items():
          if getattr(self, protKey).get():
            args += ' --{} True'.format(key)

        for key, protKey in argsEnumDic.items():
          protVal = self.getEnumText(protKey)
          args += ' --{} {}'.format(key, protVal)

        if self.fLipinski:
          lipKeys = ['LipinskiStrictFilter', 'LipinskiLenientFilter']
          args += ' --{}'.format(lipKeys[self.LipinskiType.get()])

        if self.fGhose:
          lipKeys = ['GhoseFilter', 'GhoseModifiedFilter']
          args += ' --{}'.format(lipKeys[self.GhoseType.get()])

        return args

    def getGenDirs(self, pocketDir, genMin):
      gDirs = []
      runDir = os.path.join(pocketDir, 'Run_0')
      for gDir in os.listdir(runDir):
        if 'generation_' in gDir:
          gen = int(gDir.split('generation_')[-1])
          if gen >= genMin:
            gDirs.append(os.path.join(runDir, gDir, 'PDBs'))
      return gDirs

    def parseModelEnergy(self, modelStr):
        for line in modelStr.split('\n'):
            if line.startswith('REMARK VINA RESULT:'):
              return float(line.split()[3])

    def performEnergyFilter(self, dockDic, maxEnergy):
      '''Remove those docked models with a higher energy than the specified threshold'''
      for model in dockDic:
        energy = self.parseModelEnergy(dockDic[model])
        if energy > maxEnergy:
          del dockDic[model]
      return dockDic

    def compressResults(self, genDir):
      args = '--compress_or_decompress compress --input_folder_or_file {}'.format(os.path.abspath(genDir))
      Plugin.runScript(self, 'accessory_scripts/file_concatenate_and_compression.py', args, env='AutoGrow4',
                                cwd=self._getExtraPath(),
                                scriptDir=Plugin.getProgramHome(AGROW_DIC, path='autogrow4'))

    def getDockedDic(self, genDir):
      '''Returns a dictionary with all the docked poses in a generation dir, cathegorized by docking file
      and model'''
      dockedDic = {}
      for vinaFile in os.listdir(genDir):
        if 'pdbqt.vina' in vinaFile:
          vinaFile = os.path.join(genDir, vinaFile)
          dockedDic.update(self.parsePDBQTVina(vinaFile))
      return dockedDic

    def parsePDBQTVina(self, vinaFile):
      '''Return a dictionary with all the models inside a pdbqt.vina file {model nº: pdbqtStr}'''
      modelDic = {}
      curModel, curStr = '', ''
      with open(vinaFile) as f:
        for line in f:
          if line.startswith('MODEL'):
            if curModel:
              modelBaseName = os.path.basename(vinaFile).replace('.pdbqt.vina', '')
              modelDic['{}_{}'.format(modelBaseName, curModel.split()[-1])] = curStr
            curModel, curStr = line.strip(), ''
          else:
            curStr += line
      return modelDic

    def convertInputSMIFile(self):
        '''Convert source compounds to SMI and store them into a file'''
        oDir = self.getInputLigandsPath()
        if not os.path.exists(oDir):
            os.mkdir(oDir)
        if self.ligandsSource.get() == 1:
          allMols = self.getLigandsFileNames()
          runOpenBabel(self, '{} -O {}'.format(' '.join(allMols), self.getInputSMIFile()),
                       popen=True)
        else:
          inSMIFile = Plugin.getProgramHome(AGROW_DIC, path='autogrow4/source_compounds/{}'.
                                                    format(zincDic[self.getEnumText('inputAGrowMols')]))
          shutil.copy(inSMIFile, self.getInputSMIFile())

    def getOriginalReceptorFile(self):
        if self.fromReceptor.get() == 0:
            return self.inputAtomStruct.get().getFileName()
        else:
            return self.inputStructROIs.get().getProteinFile()

    def getReceptorName(self):
        fnReceptor = self.getOriginalReceptorFile()
        return getBaseFileName(fnReceptor)

    def getReceptorDir(self):
        fnReceptor = self.getOriginalReceptorFile()
        return os.path.dirname(fnReceptor)

    def getReceptorPDBQT(self):
        return os.path.abspath(self._getExtraPath('{}.pdbqt'.format(self.getReceptorName())))

    def getReceptorPDB(self):
        return os.path.abspath(self._getExtraPath('{}.pdb'.format(self.getReceptorName())))

    def convertReceptor2PDB(self, proteinFile):
        inName, inExt = os.path.splitext(os.path.basename(proteinFile))
        oFile = self.getReceptorPDB()
        if not os.path.exists(oFile):
          args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(proteinFile), oFile)
          runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())
        return oFile

    def getConversionArgs(self):
        mglPath = Plugin.getDefPath({'name': 'mgltools', 'version': '1.5.6'})
        if os.path.exists(mglPath):
            args = ' --conversion_choice MGLToolsConversion --mgltools_directory {}'.format(mglPath)
        else:
            obabelPath = Plugin.getCondaEnvPath(env='plip', path='bin/obabel')
            args = ' --conversion_choice ObabelConversion --obabel_path {}'.format(obabelPath)
        return args

    def getOutputPocketDir(self, pocket=None):
        if pocket == None:
          outDir = os.path.abspath(self._getExtraPath('pocket_1'))
        else:
          outDir = os.path.abspath(self._getExtraPath('pocket_{}'.format(pocket.getObjId())))
        return outDir

    def getGridId(self, outDir):
        return outDir.split('_')[-1]

    def getPocketDirs(self):
        dirs = []
        for file in os.listdir(self._getExtraPath()):
            d = self._getExtraPath(file)
            if os.path.isdir(d) and 'pocket' in file:
                dirs.append(d)
        dirs.sort()
        return dirs

    def _validate(self):
      vals = []
      if self.getEnumText('scoreChoice') != 'VINA':
        mglPath = Plugin.getDefPath({'name': 'mgltools', 'version': '1.5.6'})
        if self.getEnumText('dockChoice') != 'VinaDocking' or not os.path.exists(mglPath):
          vals.append('In order to use NN1 or NN2 scoring functions, Vina docking must be selected and conversions '
                      'must be performed using MGLTools 1.5.6')
      return vals

    def _warnings(self):
      warns = []
      if self.fromReceptor.get() == 0 or len(self.inputStructROIs.get()) > 2:
          warns.append('AutoGrow4 performs several rounds of docking and genetic operations through a number of '
                       'generations, leading to a extensive computation. Using a complete protein or too many '
                       'structural regions of interest may produce very long waiting times.')
      return warns
