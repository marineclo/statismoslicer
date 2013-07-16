
// Select Parameters Widgets includes
#include "qSlicerLoadableSSMBuildingSelectParametersWidget.h"
#include "ui_qSlicerLoadableSSMBuildingSelectParametersWidget.h"

// Qt includes
#include <QFileDialog>

// Statismo includes
/*#include "statismo_ITK/itkStatisticalModel.h"
#include "itkMesh.h"
#include "itkMeshRepresenter.h"*/

#include "statismo/StatisticalModel.h"
#include "Representers/VTK/vtkPolyDataRepresenter.h"
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include "vtkSlicerModelsLogic.h"
#include "vtkSlicerLoadableSSMBuildingLogic.h"

#include <vtkNew.h>
// MRML includes
#include "vtkMRMLScene.h"
#include "vtkMRMLModelNode.h"
#include "vtkMRMLModelDisplayNode.h"

// Slicer includes
#include <qSlicerCoreApplication.h>
#include <qSlicerModuleManager.h>
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
  /*const unsigned Dimensions = 3;
  typedef itk::Mesh<float, Dimensions  > MeshType;
  typedef itk::MeshRepresenter<float, Dimensions> RepresenterType;
  typedef itk::StatisticalModel<RepresenterType> StatisticalModelType;
  typedef vnl_vector<statismo::ScalarType> VectorType;

  Q_D(qSlicerLoadableSSMBuildingSelectParametersWidget);
  // Get the model name selected by the user
  QString modelQString = d->modelNamePath->text();
  std::string modelString = modelQString.toStdString(); 
  // Load the model
  StatisticalModelType::Pointer statModel = StatisticalModelType::New();
  statModel->Load(modelString.c_str());
	int nbPCA = statModel->GetNumberOfPrincipalComponents();
  VectorType coefficients(nbPCA,0.0); // set the vector to 0

  // Compute the sample
  int nb = static_cast<int>(d->pcSlider->value());
  coefficients(nb)=d->stdSlider->value();
  //coefficients(d->pcSlider->value)=d->stdSlider->value;
  MeshType::Pointer sample = statModel->DrawSample(coefficients);*/

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
  auto_ptr<StatisticalModelType> model(StatisticalModelType::Load(modelString.c_str()));
  VectorType coefficients = VectorType::Zero(model->GetNumberOfPrincipalComponents());
  int pc = static_cast<int>(d->pcSlider->value());
  coefficients(pc) = static_cast<int>(d->stdSlider->value());
  vtkPolyData* samplePC = model->DrawSample(coefficients);
  
  std::cout << "modelString: "<<modelString<<std::endl;
  // Add polydata sample to the scene
  //vtkMRMLScene* mrmlScene = vtkMRMLScene::New();
/*
  vtkMRMLModelNode* sampleNode = vtkMRMLModelNode::New();
  sampleNode->SetScene(mrmlScene);
  sampleNode->SetName("Sample");
  sampleNode->SetAndObservePolyData(samplePC);
  
  vtkMRMLModelDisplayNode* modelDisplayNode = vtkMRMLModelDisplayNode::New();
  modelDisplayNode->SetColor(0,1,0); // green
  modelDisplayNode->SetScene(mrmlScene);
  mrmlScene->AddNode(modelDisplayNode);
  sampleNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());
  
  modelDisplayNode->SetInputPolyData(sampleNode->GetPolyData());
  mrmlScene->AddNode(sampleNode);*/

  /*vtkNew<vtkMRMLModelDisplayNode> display;
  mrmlScene->AddNode(display.GetPointer());

  vtkNew<vtkMRMLModelNode> modelDisplay;
  modelDisplay->SetAndObservePolyData(samplePC);
  modelDisplay->SetAndObserveDisplayNodeID(display->GetID());
  mrmlScene->AddNode(modelDisplay.GetPointer());*/

  // Save polydata
  vtkSlicerLoadableSSMBuildingLogic* moduleLogic = vtkSlicerLoadableSSMBuildingLogic::New();
  moduleLogic->DisplaySampleModel(samplePC);
  //qSlicerAbstractCoreModule* modelsModule = qSlicerCoreApplication::application()->moduleManager()->module("Models");
  //vtkSlicerModelsLogic* modelsLogic = vtkSlicerModelsLogic::SafeDownCast(modelsModule->logic());
  /*vtkSlicerModelsLogic* modelsLogic = vtkSlicerModelsLogic::New();
  vtkMRMLModelNode* modelNode = modelsLogic->AddModel(samplePC);
  if (modelNode == NULL){
    std::cout<<"empty"<<std::endl;
  }*/
  /*vtkPolyDataWriter *pdWriter = vtkPolyDataWriter::New();
  pdWriter->SetFileName("/home/marine/Desktop/m.vtk");
  pdWriter->SetInput(samplePC);
  pdWriter->Write();*/

    //pdWriter->Delete();


}

//-----------------------------------------------------------------------------
qSlicerLoadableSSMBuildingSelectParametersWidget
::~qSlicerLoadableSSMBuildingSelectParametersWidget()
{
}
