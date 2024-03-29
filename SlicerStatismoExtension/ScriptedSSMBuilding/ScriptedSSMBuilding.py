import os
import unittest
from __main__ import vtk, qt, ctk, slicer
import numpy as np
import statismo

#
# ScriptedSSMBuilding
#

class ScriptedSSMBuilding:
  def __init__(self, parent):
    parent.title = "ScriptedSSMBuilding" 
    parent.categories = ["Examples"]
    parent.dependencies = []
    parent.contributors = ["Marine Clogenson (EPFL)"] 
    parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    """
    parent.acknowledgementText = """
    This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc. and Steve Pieper, Isomics, Inc.  and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.
    self.parent = parent


#
# qScriptedSSMBuildingWidget
#

class ScriptedSSMBuildingWidget:
  def __init__(self, parent = None):
    if not parent:
      self.parent = slicer.qMRMLWidget()
      self.parent.setLayout(qt.QVBoxLayout())
      self.parent.setMRMLScene(slicer.mrmlScene)
    else:
      self.parent = parent
    self.layout = self.parent.layout()
    if not parent:
      self.setup()
      self.parent.show()

  def setup(self):

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    #  Import model
    #
    self.inputModel = qt.QFileDialog.getOpenFileName()
    parametersFormLayout.addRow("Input Model: ", qt.QLabel(self.inputModel))

    # Principal Component selection scroller
    self.pcSlider = ctk.ctkSliderWidget()
    self.pcSlider.decimals = 0
    self.pcSlider.minimum = 0
    self.pcSlider.maximum = 10
    parametersFormLayout.addRow("PC", self.pcSlider)
    
    # Standard variation selection scroller
    self.stdSlider = ctk.ctkSliderWidget()
    self.stdSlider.decimals = 0
    self.stdSlider.minimum = -10
    self.stdSlider.maximum = 10
    parametersFormLayout.addRow("std", self.stdSlider)

    # make connections
    self.pcSlider.connect('valueChanged(double)', self.onSelect)
    self.stdSlider.connect('valueChanged(double)', self.onSelect)

    # make an instance of the logic for use by the slots
    self.logic = ScriptedSSMBuildingLogic()

    # Add vertical spacer
    self.layout.addStretch(1)

  def onSelect(self):
    print self.inputModel
    print self.pcSlider.value
    print self.stdSlider.value
    self.logic.selectSample(self.inputModel, self.pcSlider.value, self.stdSlider.value)
    
  def cleanup(self):
    pass

#
# ScriptedSSMBuildingLogic
#

class ScriptedSSMBuildingLogic:
  """This class should implement all the actual 
  computation done by your module.  The interface 
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget
  """
  def __init__(self):
    pass

  def selectSample(self,model, pcIndex,stdIndex):
    #slicer.util.loadModel(model)
    # load the model
    stat_model = statismo.StatisticalModel_vtkPD.Load(model)
    nbPC = stat_model.GetNumberOfPrincipalComponents()
    coefficients = np.zeros(nbPC)
    #Compute the Sample
    coefficients[pcIndex] = stdIndex
    samplePC = stat_model.DrawSample(coefficients)
    
    # Add polydata to the scene 
    # Need to be tested !!!
    modelNode = slicer.vtkMRMLNode()
    modelNode.SetPolyData(samplePC)
    slicer.mrmlScene.AddNode(modelNode)
    
    print model
    print pcIndex
    print stdIndex

