# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from pwem.protocols import ProtImportPdb, ProtSetFilter

from pwchem.tests import TestExtractLigand, TestImportBase, TestDefineStructROIs

from durrantlab.protocols import *

defROIsStr = '''1) Ligand: {"molName": "HEM"}'''

class TestGypSumDL(TestImportBase):
    @classmethod
    def _runGypsumDLPrep(cls, inProt):
        protGyp = cls.newProtocol(
            ProtChemGypsumDL,
        )
        protGyp.inputSmallMolecules.set(inProt)
        protGyp.inputSmallMolecules.setExtended('outputSmallMolecules')

        cls.proj.launchProtocol(protGyp, wait=True)
        return protGyp

    def test(self):
        protGyp = self._runGypsumDLPrep(inProt=self.protImportSmallMols)
        self.assertIsNotNone(getattr(protGyp, 'outputSmallMolecules', None))

class TestAutoGrow(TestDefineStructROIs):
    @classmethod
    def _runAutoGrow(cls, inProt):
        protGyp = cls.newProtocol(
            ProtChemAutoGrow4,
            fromReceptor=1, ligandsSource=0,
            nGens=2, nCrossovers=3, nMutants=3, nElitism=3, topMolsSeed=3,
        )
        protGyp.inputStructROIs.set(inProt)
        protGyp.inputStructROIs.setExtended('outputStructROIs')

        cls.proj.launchProtocol(protGyp, wait=True)
        return protGyp

    def _runFilterSites(self, inProt):
        protFilter = self.newProtocol(
            ProtSetFilter,
            operation=ProtSetFilter.CHOICE_RANKED,
            threshold=1,
            rankingField='_score')
        protFilter.inputSet.set(inProt)
        protFilter.inputSet.setExtended('outputStructROIs')

        self.launchProtocol(protFilter)
        return protFilter

    def test(self):
        pDef = self._runDefStructROIs(defROIsStr)
        self._waitOutput(pDef, 'outputStructROIs', sleepTime=10)
        pROI = self._runFilterSites(pDef)

        pAGrow = self._runAutoGrow(pROI)
        self.assertIsNotNone(getattr(pAGrow, 'outputSmallMolecules', None))

class TestDeepFrag(TestExtractLigand):
    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0, pdbId='4erf')
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runDeepFrag(cls, inProt):
        protGyp = cls.newProtocol(
            ProtChemDeepFrag,
            inputLigand='SmallMolecule (g1_4erf_0R3-1_1 molecule)', topK=10,
            ligandList='{"Ligand": "SmallMolecule (g1_4erf_0R3-1_1 molecule)", "ConnectionPoint": "C23"}'
        )
        protGyp.inputSmallMolecules.set(inProt)
        protGyp.inputSmallMolecules.setExtended('outputSmallMolecules')

        cls.proj.launchProtocol(protGyp, wait=True)
        return protGyp

    def test(self):
        protExtract = self._runExtractLigand(self.protImportPDB)
        self._waitOutput(protExtract, 'outputSmallMolecules', sleepTime=5)

        protDeepFrag = self._runDeepFrag(inProt=protExtract)
        self.assertIsNotNone(protDeepFrag.outputSmallMolecules,
                             "There was a problem with DeepFrag")
