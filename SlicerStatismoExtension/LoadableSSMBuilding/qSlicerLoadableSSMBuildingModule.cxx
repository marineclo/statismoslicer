
// Qt includes
#include <QtPlugin>

// LoadableSSMBuilding Logic includes
#include <vtkSlicerLoadableSSMBuildingLogic.h>

// LoadableSSMBuilding includes
#include "qSlicerLoadableSSMBuildingModule.h"
#include "qSlicerLoadableSSMBuildingModuleWidget.h"

//-----------------------------------------------------------------------------
Q_EXPORT_PLUGIN2(qSlicerLoadableSSMBuildingModule, qSlicerLoadableSSMBuildingModule);

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerLoadableSSMBuildingModulePrivate
{
public:
  qSlicerLoadableSSMBuildingModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerLoadableSSMBuildingModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerLoadableSSMBuildingModulePrivate
::qSlicerLoadableSSMBuildingModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerLoadableSSMBuildingModule methods

//-----------------------------------------------------------------------------
qSlicerLoadableSSMBuildingModule
::qSlicerLoadableSSMBuildingModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerLoadableSSMBuildingModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerLoadableSSMBuildingModule::~qSlicerLoadableSSMBuildingModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerLoadableSSMBuildingModule::helpText()const
{
  return "This is a loadable module for Atlas building";
}

//-----------------------------------------------------------------------------
QString qSlicerLoadableSSMBuildingModule::acknowledgementText()const
{
  //return "This work was was partially funded by NIH grant 3P41RR013218-12S1";
  return "";
}

//-----------------------------------------------------------------------------
QStringList qSlicerLoadableSSMBuildingModule::contributors()const
{
  QStringList moduleContributors;
  moduleContributors << QString("Marine Clogenson (EPFL)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerLoadableSSMBuildingModule::icon()const
{
  //return QIcon(":/Icons/LoadableSSMBuilding.png");
  return QIcon("");
}

//-----------------------------------------------------------------------------
QStringList qSlicerLoadableSSMBuildingModule::categories() const
{
  return QStringList() << "Examples";
}

//-----------------------------------------------------------------------------
QStringList qSlicerLoadableSSMBuildingModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation * qSlicerLoadableSSMBuildingModule
::createWidgetRepresentation()
{
  return new qSlicerLoadableSSMBuildingModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerLoadableSSMBuildingModule::createLogic()
{
  return vtkSlicerLoadableSSMBuildingLogic::New();
}
