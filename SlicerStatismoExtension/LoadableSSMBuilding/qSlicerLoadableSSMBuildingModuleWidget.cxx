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
#include <QDebug>
#include <QFileDialog>

// SlicerQt includes
#include "qSlicerLoadableSSMBuildingModuleWidget.h"
#include "ui_qSlicerLoadableSSMBuildingModuleWidget.h"

// Statismo includes
#include "statismo/StatisticalModel.h"
#include "Representers/VTK/vtkPolyDataRepresenter.h"
#include <vtkPolyData.h>
//#include <vtkPolyDataWriter.h>

// LoadableSSMBuilding Logic includes
#include "vtkSlicerLoadableSSMBuildingLogic.h"

#include "vtkMRMLModelNode.h"
#include "vtkMRMLModelDisplayNode.h"
#include "vtkMRMLScene.h"
#include "vtkMRMLLayoutNode.h"
#include "vtkMRMLDoubleArrayNode.h"
#include <vtkDoubleArray.h>
#include "vtkMRMLChartViewNode.h"
#include "vtkMRMLChartNode.h"


using namespace statismo;
typedef vtkPolyDataRepresenter RepresenterType;
typedef StatisticalModel<RepresenterType> StatisticalModelType;

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerLoadableSSMBuildingModuleWidgetPrivate: public Ui_qSlicerLoadableSSMBuildingModuleWidget
{
public:
  qSlicerLoadableSSMBuildingModuleWidgetPrivate();

  StatisticalModelType* statModel;
  void setStatModel(StatisticalModelType* newStatModel);
  StatisticalModelType* getStatModel();
};

//-----------------------------------------------------------------------------
// qSlicerLoadableSSMBuildingModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerLoadableSSMBuildingModuleWidgetPrivate::qSlicerLoadableSSMBuildingModuleWidgetPrivate()
{
}

void qSlicerLoadableSSMBuildingModuleWidgetPrivate::setStatModel(StatisticalModelType* newStatModel)
{
  this->statModel = newStatModel;
}

StatisticalModelType* qSlicerLoadableSSMBuildingModuleWidgetPrivate::getStatModel()
{
  StatisticalModelType* newStatModel;
  return newStatModel = this->statModel;
}
//-----------------------------------------------------------------------------
// qSlicerLoadableSSMBuildingModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerLoadableSSMBuildingModuleWidget::qSlicerLoadableSSMBuildingModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerLoadableSSMBuildingModuleWidgetPrivate )
{
}

//-----------------------------------------------------------------------------
qSlicerLoadableSSMBuildingModuleWidget::~qSlicerLoadableSSMBuildingModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::setup()
{
  Q_D(qSlicerLoadableSSMBuildingModuleWidget);
  d->setupUi(this);
  
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::onSelectInputModel()
{
  Q_D(qSlicerLoadableSSMBuildingModuleWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select input model", QString());
  d->modelNamePath->setText(inputFile);
  // Load the model
  std::string modelString = inputFile.toStdString(); 
  d->setStatModel(d->statModel->Load(modelString.c_str()));
  // Display Eigen spectrum
  displayEigenSpectrum();
  // Set the number of components for the pcSlider
  d->pcSlider->setMaximum(d->statModel->GetNumberOfPrincipalComponents());

  //Calculate mean
  vtkPolyData* meanModel = d->statModel->DrawMean();
  // Add mean model to the scene
  vtkSlicerLoadableSSMBuildingLogic* moduleLogic = vtkSlicerLoadableSSMBuildingLogic::New();
  moduleLogic->DisplaySampleModel(meanModel, this->mrmlScene());
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::onSelect()
{
  /**** VTK Model ****/
  using std::auto_ptr;
  Q_D(qSlicerLoadableSSMBuildingModuleWidget);
  // Get the model name selected by the user
  d->statModel = d->getStatModel();
  int nbPrincipalComponent = d->statModel->GetNumberOfPrincipalComponents();
  VectorType coefficients = VectorType::Zero(nbPrincipalComponent);
  int pc = static_cast<int>(d->pcSlider->value())-1; // -1 bacause user can choose between 1 and max
  coefficients(pc) = d->stdSlider->value();
  vtkPolyData* samplePC = d->statModel->DrawSample(coefficients);
  
  // Add polydata sample to the scene
  vtkSlicerLoadableSSMBuildingLogic* moduleLogic = vtkSlicerLoadableSSMBuildingLogic::New();
  moduleLogic->DisplaySampleModel(samplePC, this->mrmlScene());
  /*vtkMRMLScene* mrmlScene = this->mrmlScene();
  vtkMRMLModelNode* sampleNode = vtkMRMLModelNode::New();
  sampleNode->SetScene(mrmlScene);
  sampleNode->SetName("Sample");
  sampleNode->SetAndObservePolyData(samplePC);
  
  vtkMRMLModelDisplayNode* modelDisplayNode = vtkMRMLModelDisplayNode::New();
  //modelDisplayNode->SetColor(0,1,0); // green
  modelDisplayNode->SetScene(mrmlScene);
  mrmlScene->AddNode(modelDisplayNode);
  sampleNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());
  
  modelDisplayNode->SetInputPolyData(sampleNode->GetPolyData());
  mrmlScene->AddNode(sampleNode);*/
  
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::setMRMLScene(vtkMRMLScene* mrmlScene)
{
  this->qSlicerAbstractModuleWidget::setMRMLScene(mrmlScene);
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::displayEigenSpectrum()
{
  Q_D(qSlicerLoadableSSMBuildingModuleWidget);
  // Get the model name selected by the user
  d->statModel = d->getStatModel();

  // Switch to a layout (24) that contains a Chart View to initiate the construction of the widget and Chart View Node 
  vtkMRMLScene* mrmlScene = this->mrmlScene();
  mrmlScene->InitTraversal();
  vtkMRMLLayoutNode* sceneLayoutNode = vtkMRMLLayoutNode::SafeDownCast(mrmlScene->GetNextNodeByClass("vtkMRMLLayoutNode") );
  sceneLayoutNode->SetViewArrangement(24);

  // Get the Chart View Node
  vtkMRMLChartViewNode* chartViewNode = vtkMRMLChartViewNode::SafeDownCast(mrmlScene->GetNextNodeByClass("vtkMRMLChartViewNode") );

  // Create an Array Node and add the eigen value
  VectorType eigenvalue = d->statModel->GetPCAVarianceVector();
  vtkMRMLDoubleArrayNode* doubleArrayNode = vtkMRMLDoubleArrayNode::SafeDownCast(mrmlScene->AddNode(vtkMRMLDoubleArrayNode::New()));
  vtkDoubleArray* a = doubleArrayNode->GetArray();
  int nbPrincipalComponent = d->statModel->GetNumberOfPrincipalComponents();
  a->SetNumberOfTuples(nbPrincipalComponent);
  for (int i = 0;i<nbPrincipalComponent;i++){
    a->SetComponent(i, 0, i);
    a->SetComponent(i, 1, eigenvalue[i]);
    a->SetComponent(i, 2, 0);
  }

  // Create a Chart Node.
  vtkMRMLChartNode* chartNode = vtkMRMLChartNode::New();
  chartNode->AddArray("EigenSpectrum", doubleArrayNode->GetID());
  mrmlScene->AddNode(chartNode);
  
  // Set properties on the Chart
  chartNode->SetProperty("default", "title", "EigenSpectrum");
  chartNode->SetProperty("default", "xAxisLabel", "Principal components");
  chartNode->SetProperty("default", "yAxisLabel", "Eigen value");
  chartNode->SetProperty("default", "type", "Line");
  chartNode->SetProperty("default", "showMarkers", "on");
  
  // Tell the Chart View which Chart to display
  chartViewNode->SetChartNodeID(chartNode->GetID());
}
