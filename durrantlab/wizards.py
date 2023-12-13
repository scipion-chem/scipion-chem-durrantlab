# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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

"""
This wizard will show the structure of the pdb using a matplotlib viewer
to select the radius of the sphere that contains the protein or a desired zone.
"""

import os

from pwchem.wizards import GetRadiusProtein, SelectElementWizard, SelectLigandAtom, \
  AddElementWizard, VariableWizard
from pwchem.utils import getBaseName
from pwchem.viewers import PyMolViewer

from durrantlab.protocols import ProtChemDeepFrag, ProtChemAutoGrow4

for prot in [ProtChemAutoGrow4]:
    GetRadiusProtein().addTarget(protocol=prot,
                                 targets=['radius'],
                                 inputs=['inputAtomStruct'],
                                 outputs=['radius'])

SelectElementWizard().addTarget(protocol=ProtChemDeepFrag,
                                targets=['inputLigand'],
                                inputs=['inputSmallMolecules'],
                                outputs=['inputLigand'])

SelectLigandAtom().addTarget(protocol=ProtChemDeepFrag,
                             targets=['conPoint'],
                             inputs=['inputSmallMolecules', 'inputLigand'],
                             outputs=['conPoint'])

SelectLigandAtom().addTarget(protocol=ProtChemDeepFrag,
                             targets=['remPoint'],
                             inputs=['inputSmallMolecules', 'inputLigand'],
                             outputs=['remPoint'])


class AddDeepFragLigand(AddElementWizard):
  """Add Ligand-ConPoint-RemPoint to the list for DeepFrag execution"""
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    protocol = form.protocol

    ligName = getattr(protocol, inputParam[0]).get().strip()
    conPoint, remPoint = getattr(protocol, inputParam[1]).get().strip(), \
                         getattr(protocol, inputParam[2]).get().strip()

    if ligName != '' and conPoint != '':
      prevList = self.curePrevList(getattr(protocol, outputParam[0]).get())

      towrite = prevList + '{"Ligand": "%s", "ConnectionPoint": "%s"' % (ligName, conPoint)
      if remPoint.strip() != '':
        towrite += ', "RemovalPoint": "%s"' % remPoint
      towrite += '}\n'

      form.setVar(outputParam[0], towrite)

AddDeepFragLigand().addTarget(protocol=ProtChemDeepFrag,
                              targets=['addLigand'],
                              inputs=['inputLigand', 'conPoint', 'remPoint'],
                              outputs=['ligandList'])

class ViewInputLigandWizard(VariableWizard):
  """Visualize the chosen ligand with the correspondant labels"""
  _targets, _inputs, _outputs = [], {}, {}

  def getMol(self, inSet, molName):
    myMol = None
    for mol in inSet:
      if mol.__str__() == molName:
        myMol = mol.clone()
        break
    if myMol == None:
      print('The input ligand is not found')
      return None
    else:
      return myMol

  def writePmlFile(self, pmlFile, molFile, molName, targetFile=None):
    pmlStr = ''
    if targetFile:
      pmlStr = 'load {}\n'.format(targetFile)

    molFile = os.path.abspath(molFile)
    pmlStr += 'load {}, {}\nhide spheres, {}\nshow sticks, {}\n'. \
      format(molFile, molName, molName, molName)
    pmlStr += 'label {}, name'.format(molName)

    with open(pmlFile, 'w') as f:
      f.write(pmlStr)

  def show(self, form, *params):
    inputParam, _ = self.getInputOutput(form)
    protocol = form.protocol
    project = protocol.getProject()

    inSet, molName = getattr(protocol, inputParam[0]).get(), getattr(protocol, inputParam[1]).get()
    mol, targetFile = self.getMol(inSet, molName), inSet.getProteinFile()
    molFile = mol.getPoseFile() if mol.getPoseFile() else mol.getFileName()
    molName = getBaseName(molFile)

    pmlsDir = project.getTmpPath()
    pmlFile = os.path.join(pmlsDir, '{}.pml'.format(molName))
    self.writePmlFile(pmlFile, molFile, molName, targetFile)

    pymolV = PyMolViewer(project=project)
    view = pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))[0]
    view.show()

ViewInputLigandWizard().addTarget(protocol=ProtChemDeepFrag,
                              targets=['viewLigand'],
                              inputs=['inputSmallMolecules', 'inputLigand'],
                              outputs=[])