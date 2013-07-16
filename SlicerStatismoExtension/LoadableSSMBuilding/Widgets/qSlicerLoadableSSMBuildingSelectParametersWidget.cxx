
// Select Parameters Widgets includes
#include "qSlicerLoadableSSMBuildingSelectParametersWidget.h"
#include "ui_qSlicerLoadableSSMBuildingSelectParametersWidget.h"

// Qt includes
#include <QFileDialog>

// Statismo includes
#include "statismo/StatisticalModel.h"
#include "Representers/VTK/vtkPolyDataRepresenter.h"
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>

// LoadableSSMBuilding Logic includes
#include "vtkSlicerLoadableSSMBuildingLogic.h"


//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_LoadableSSMBuilding
class qSlicerLoadableSSMBuildingSelectParametersWidgetPrivate
  : public Ui_qSlicerLoadableSSMBuildingSelectParametersWidget
{
  Q_DECLARE_PUBLIC(qSlicerLoadableSSMBuildingSelectParametersWidget);
protected:
  qSlicerLoadableSSMBuildingSelectParametersWidget* const q_ptr;

public:
  qSlicerLoadableSSMBuildingSelectParametersWidgetPrivate(
    qSlicerLoadableSSMBuildingSelectParametersWidget& object);
  virtual void setupUi(qSlicerLoadableSSMBuildingSelectParametersWidget*);
};

// --------------------------------------------------------------------------
qSlicerLoadableSSMBuildingSelectParametersWidgetPrivate
::qSlicerLoadableSSMBuildingSelectParametersWidgetPrivate(
  qSlicerLoadableSSMBuildingSelectParametersWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingSelectParametersWidgetPrivate
::setupUi(qSlicerLoadableSSMBuildingSelectParametersWidget* widget)
{
  this->Ui_qSlicerLoadableSSMBuildingSelectParametersWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerLoadableSSMBuildingSelectParametersWidget methods

//-----------------------------------------------------------------------------
qSlicerLoadableSSMBuildingSelectParametersWidget
::qSlicerLoadableSSMBuildingSelectParametersWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerLoadableSSMBuildingSelectParametersWidgetPrivate(*this) )
{
  Q_D(qSlicerLoadableSSMBuildingSelectParametersWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingSelectParametersWidget::onSelectInputModel()
{
  Q_D(qSlicerLoadableSSMBuildingSelectParametersWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select input model", QString());
  d->modelNamePath->setText(inputFile);

}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingSelectParametersWidget::onSelect()
{
  /**** VTK Model ****/
  using std::auto_ptr;
  using namespace statismo;
  typedef vtkPolyDataRepresenter RepresenterType;
  typedef StatisticalModel<RepresenterType> StatisticalModelType;
  Q_D(qSlicerLoadableSSMBuildingSelectParametersWidget);
  // Get the model name selected by the user
  QString modelQString = d->modelNamePath->text();
  std::string modelString = modelQString.toStdString(); 
  // Load the model
  auto_ptr<StatisticalModelType> statModel(StatisticalModelType::Load(modelString.c_str()));
  VectorType coefficients = VectorType::Zero(statModel->GetNumberOfPrincipalComponents());
  int pc = static_cast<int>(d->pcSlider->value());
  coefficients(pc) = d->stdSlider->value();
  vtkPolyData* samplePC = statModel->DrawSample(coefficients);
  
  // Add polydata sample to the scene
  vtkSlicerLoadableSSMBuildingLogic* moduleLogic = vtkSlicerLoadableSSMBuildingLogic::New();
  moduleLogic->DisplaySampleModel(samplePC);
  
  // Save polydata
  /*vtkPolyDataWriter *pdWriter = vtkPolyDataWriter::New();
  pdWriter->SetFileName("/home/marine/Desktop/m.vtk");
  pdWriter->SetInput(samplePC);
  pdWriter->Write();*/

}

//-----------------------------------------------------------------------------
qSlicerLoadableSSMBuildingSelectParametersWidget
::~qSlicerLoadableSSMBuildingSelectParametersWidget()
{
}
