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

// SSMViewer Logic includes
#include <vtkSlicerSSMViewerLogic.h>

// SSMViewer includes
#include "qSlicerSSMViewerModule.h"
#include "qSlicerSSMViewerModuleWidget.h"

//-----------------------------------------------------------------------------
Q_EXPORT_PLUGIN2(qSlicerSSMViewerModule, qSlicerSSMViewerModule);

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerSSMViewerModulePrivate
{
public:
  qSlicerSSMViewerModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerSSMViewerModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerSSMViewerModulePrivate
::qSlicerSSMViewerModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerSSMViewerModule methods

//-----------------------------------------------------------------------------
qSlicerSSMViewerModule
::qSlicerSSMViewerModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerSSMViewerModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerSSMViewerModule::~qSlicerSSMViewerModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerSSMViewerModule::helpText()const
{
  return "Module for displaying a Statistical Shape Model (SSM) created with Statismo<br>";
}

//-----------------------------------------------------------------------------
QString qSlicerSSMViewerModule::acknowledgementText()const
{
  //return "This work was was partially funded by NIH grant 3P41RR013218-12S1";
  return "";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSSMViewerModule::contributors()const
{
  QStringList moduleContributors;
  moduleContributors << QString("Marine Clogenson (EPFL), Jean-Christophe Fillion-Robin (Kitware)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerSSMViewerModule::icon()const
{
  //return QIcon(":/Icons/SSMViewer.png");
  return QIcon("");
}

//-----------------------------------------------------------------------------
QStringList qSlicerSSMViewerModule::categories() const
{
  return QStringList() << "Examples";
}

//-----------------------------------------------------------------------------
QStringList qSlicerSSMViewerModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerSSMViewerModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation * qSlicerSSMViewerModule
::createWidgetRepresentation()
{
  return new qSlicerSSMViewerModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerSSMViewerModule::createLogic()
{
  return vtkSlicerSSMViewerLogic::New();
}

//-----------------------------------------------------------------------------
void qSlicerSSMViewerModule::setMRMLScene(vtkMRMLScene* mrmlScene)
{
  this->Superclass::setMRMLScene(mrmlScene);
}
