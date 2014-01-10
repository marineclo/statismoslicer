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

// LandmarkSegmentation Logic includes
#include <vtkSlicerLandmarkSegmentationLogic.h>

// LandmarkSegmentation includes
#include "qSlicerLandmarkSegmentationModule.h"
#include "qSlicerLandmarkSegmentationModuleWidget.h"

//-----------------------------------------------------------------------------
Q_EXPORT_PLUGIN2(qSlicerLandmarkSegmentationModule, qSlicerLandmarkSegmentationModule);

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerLandmarkSegmentationModulePrivate
{
public:
  qSlicerLandmarkSegmentationModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerLandmarkSegmentationModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerLandmarkSegmentationModulePrivate
::qSlicerLandmarkSegmentationModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerLandmarkSegmentationModule methods

//-----------------------------------------------------------------------------
qSlicerLandmarkSegmentationModule
::qSlicerLandmarkSegmentationModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerLandmarkSegmentationModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerLandmarkSegmentationModule::~qSlicerLandmarkSegmentationModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerLandmarkSegmentationModule::helpText()const
{
  return "This is a loadable module bundled in an extension";
}

//-----------------------------------------------------------------------------
QString qSlicerLandmarkSegmentationModule::acknowledgementText()const
{
  return "This work was was partially funded by NIH grant 3P41RR013218-12S1";
}

//-----------------------------------------------------------------------------
QStringList qSlicerLandmarkSegmentationModule::contributors()const
{
  QStringList moduleContributors;
  moduleContributors << QString("Jean-Christophe Fillion-Robin (Kitware)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerLandmarkSegmentationModule::icon()const
{
  return QIcon(":/Icons/LandmarkSegmentation.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerLandmarkSegmentationModule::categories() const
{
  return QStringList() << "Examples";
}

//-----------------------------------------------------------------------------
QStringList qSlicerLandmarkSegmentationModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerLandmarkSegmentationModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation * qSlicerLandmarkSegmentationModule
::createWidgetRepresentation()
{
  return new qSlicerLandmarkSegmentationModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerLandmarkSegmentationModule::createLogic()
{
  return vtkSlicerLandmarkSegmentationLogic::New();
}
