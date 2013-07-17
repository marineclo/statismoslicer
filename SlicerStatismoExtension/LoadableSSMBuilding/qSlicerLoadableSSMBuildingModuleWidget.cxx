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

  //StatisticalModelType* statModel;
};

//-----------------------------------------------------------------------------
// qSlicerLoadableSSMBuildingModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerLoadableSSMBuildingModuleWidgetPrivate::qSlicerLoadableSSMBuildingModuleWidgetPrivate()
{
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
  //std::string modelString = inputFile.toStdString(); 
  //d->statModel->Load(modelString.c_str());
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::onSelect()
{
  /**** VTK Model ****/
  using std::auto_ptr;
  Q_D(qSlicerLoadableSSMBuildingModuleWidget);
  // Get the model name selected by the user
  QString modelQString = d->modelNamePath->text();
  std::string modelString = modelQString.toStdString(); 
  // Load the model
  auto_ptr<StatisticalModelType> statModel(StatisticalModelType::Load(modelString.c_str()));
  int nbPrincipalComponent = statModel->GetNumberOfPrincipalComponents();
  VectorType coefficients = VectorType::Zero(nbPrincipalComponent);
  int pc = static_cast<int>(d->pcSlider->value());
  coefficients(pc) = d->stdSlider->value();
  vtkPolyData* samplePC = statModel->DrawSample(coefficients);
  
  //std::cout<<"output "<<this->mrmlScene()<<std::endl;
  // Add polydata sample to the scene
  //vtkSlicerLoadableSSMBuildingLogic* moduleLogic = vtkSlicerLoadableSSMBuildingLogic::New();
  //moduleLogic->DisplaySampleModel(samplePC);
  vtkMRMLScene* mrmlScene = this->mrmlScene();
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
  mrmlScene->AddNode(sampleNode);
  
  // Switch to a layout (24) that contains a Chart View to initiate the construction of the widget and Chart View Node 
  mrmlScene->InitTraversal();
  vtkMRMLLayoutNode* sceneLayoutNode = vtkMRMLLayoutNode::SafeDownCast(mrmlScene->GetNextNodeByClass("vtkMRMLLayoutNode") );
  sceneLayoutNode->SetViewArrangement(24);

  // Get the Chart View Node
  vtkMRMLChartViewNode* chartViewNode = vtkMRMLChartViewNode::SafeDownCast(mrmlScene->GetNextNodeByClass("vtkMRMLChartViewNode") );
  if (chartViewNode == NULL)
  {
    std::cout<<"output chartViewNode "<<std::endl;
  }

  // Create an Array Node and add the eigen value
  VectorType eigenvalue = statModel->GetPCAVarianceVector();
  vtkMRMLDoubleArrayNode* doubleArrayNode = vtkMRMLDoubleArrayNode::New();
  vtkDoubleArray* a = vtkDoubleArray::New();
  a->SetNumberOfTuples(nbPrincipalComponent);
  std::cout<<"output "<<nbPrincipalComponent<<std::endl;
  for (int i = 0;i<nbPrincipalComponent;i++){
    //std::cout<<"output "<<eigenvalue[i]<<std::endl;
    a->SetComponent(i, 0, i);
    a->SetComponent(i, 1, eigenvalue[i]);
    a->SetComponent(i, 2, 0);
  }
  doubleArrayNode->SetArray(a);
  mrmlScene->AddNode(doubleArrayNode);

  // Create a Chart Node.
  vtkMRMLChartNode* chartNode = vtkMRMLChartNode::New();
  chartNode->AddArray("EigenSpectrum", doubleArrayNode->GetID());
  mrmlScene->AddNode(chartNode);
  
  // Set a few properties on the Chart. The first argument is a string identifying which Array to assign the property. 
  chartNode->SetProperty("default", "title", "A simple chart with 2 curves");
  chartNode->SetProperty("default", "xAxisLabel", "Something in x");
  chartNode->SetProperty("default", "yAxisLabel", "Something in y");
  

  // Tell the Chart View which Chart to display
  //vtkMRMLChartViewNode* chartViewNode = vtkMRMLChartViewNode::New();
  chartViewNode->SetChartNodeID(chartNode->GetID());
  //mrmlScene->AddNode(chartViewNode);

  
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::setMRMLScene(vtkMRMLScene* mrmlScene)
{
  this->qSlicerAbstractModuleWidget::setMRMLScene(mrmlScene);
}

