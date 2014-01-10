/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// Qt includes
#include <QtPlugin>

// DisplaySSM Logic includes
#include <vtkSlicerDisplaySSMLogic.h>

// DisplaySSM includes
#include "qSlicerDisplaySSMModule.h"
#include "qSlicerDisplaySSMModuleWidget.h"

//-----------------------------------------------------------------------------
Q_EXPORT_PLUGIN2(qSlicerDisplaySSMModule, qSlicerDisplaySSMModule);

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerDisplaySSMModulePrivate
{
public:
  qSlicerDisplaySSMModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerDisplaySSMModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerDisplaySSMModulePrivate
::qSlicerDisplaySSMModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerDisplaySSMModule methods

//-----------------------------------------------------------------------------
qSlicerDisplaySSMModule
::qSlicerDisplaySSMModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerDisplaySSMModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerDisplaySSMModule::~qSlicerDisplaySSMModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerDisplaySSMModule::helpText()const
{
  return "Module for displaying a Statistical Shape Model (SSM) created with Statismo<br>";
}

//-----------------------------------------------------------------------------
QString qSlicerDisplaySSMModule::acknowledgementText()const
{
  //return "This work was was partially funded by NIH grant 3P41RR013218-12S1";
  return "";
}

//-----------------------------------------------------------------------------
QStringList qSlicerDisplaySSMModule::contributors()const
{
  QStringList moduleContributors;
  moduleContributors << QString("Marine Clogenson (EPFL), Jean-Christophe Fillion-Robin (Kitware)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerDisplaySSMModule::icon()const
{
  //return QIcon(":/Icons/DisplaySSM.png");
  return QIcon("");
}

//-----------------------------------------------------------------------------
QStringList qSlicerDisplaySSMModule::categories() const
{
  return QStringList() << "Examples";
}

//-----------------------------------------------------------------------------
QStringList qSlicerDisplaySSMModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerDisplaySSMModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation * qSlicerDisplaySSMModule
::createWidgetRepresentation()
{
  return new qSlicerDisplaySSMModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerDisplaySSMModule::createLogic()
{
  return vtkSlicerDisplaySSMLogic::New();
}

//-----------------------------------------------------------------------------
void qSlicerDisplaySSMModule::setMRMLScene(vtkMRMLScene* mrmlScene)
{
  this->Superclass::setMRMLScene(mrmlScene);
}
