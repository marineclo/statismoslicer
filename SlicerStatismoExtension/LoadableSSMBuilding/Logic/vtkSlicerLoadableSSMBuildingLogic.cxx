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

// LoadableSSMBuilding Logic includes
#include "vtkSlicerLoadableSSMBuildingLogic.h"

// MRML includes

// VTK includes
#include <vtkNew.h>
#include "vtkMRMLModelNode.h"
#include "vtkMRMLModelDisplayNode.h"

#include "vtkMRMLLayoutNode.h"
#include "vtkMRMLDoubleArrayNode.h"
#include <vtkDoubleArray.h>
#include "vtkMRMLChartViewNode.h"
#include "vtkMRMLChartNode.h"

// STD includes
#include <cassert>

#include "vtkSlicerModelsLogic.h"

#include "qSlicerAbstractCoreModule.h"
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerLoadableSSMBuildingLogic);

//----------------------------------------------------------------------------
vtkSlicerLoadableSSMBuildingLogic::vtkSlicerLoadableSSMBuildingLogic()
{
}

//----------------------------------------------------------------------------
vtkSlicerLoadableSSMBuildingLogic::~vtkSlicerLoadableSSMBuildingLogic()
{
}

//----------------------------------------------------------------------------
void vtkSlicerLoadableSSMBuildingLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerLoadableSSMBuildingLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerLoadableSSMBuildingLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerLoadableSSMBuildingLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerLoadableSSMBuildingLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerLoadableSSMBuildingLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerLoadableSSMBuildingLogic::SetModelsLogic(vtkSlicerModelsLogic* modelsLogic)
{
}

//---------------------------------------------------------------------------
void vtkSlicerLoadableSSMBuildingLogic::DisplaySampleModel(vtkPolyData* polydata, vtkMRMLScene * mrmlScene, StatisticalModelType* statModel)
{
  this->SetMRMLScene(mrmlScene);
  vtkSlicerModelsLogic* modelsLogic = vtkSlicerModelsLogic::New();
  modelsLogic->SetMRMLScene(mrmlScene);
  vtkMRMLModelNode* modelNode = modelsLogic->AddModel(polydata);
  if (modelNode == NULL){
    std::cout<<"empty"<<std::endl;
  }
  mrmlScene = this->GetMRMLScene();
  std::cout<<"output "<<this->GetMRMLScene()<<std::endl;

// Switch to a layout (24) that contains a Chart View to initiate the construction of the widget and Chart View Node 
  /*vtkSmartPointer<vtkCollection> layoutNodes = vtkSmartPointer<vtkCollection>::Take( mrmlScene->GetNodesByClass("vtkMRMLLayoutNode") );
  layoutNodes->InitTraversal();
  vtkObject* layoutNodeVtkObject = layoutNodes->GetNextItemAsObject();
  vtkMRMLLayoutNode* layoutNode = vtkMRMLLayoutNode::SafeDownCast(layoutNodeVtkObject);
  layoutNode->SetViewArrangement(24);*/
  
  mrmlScene->InitTraversal();
  vtkMRMLLayoutNode* sceneLayoutNode = vtkMRMLLayoutNode::SafeDownCast(mrmlScene->GetNextNodeByClass("vtkMRMLLayoutNode") );
  sceneLayoutNode->SetViewArrangement(24);

  /*vtkSmartPointer<vtkCollection> chartViewNodes = vtkSmartPointer<vtkCollection>::Take( mrmlScene->GetNodesByClass("vtkMRMLChartViewNode") );
  chartViewNodes->InitTraversal();
  vtkObject* chartViewNodeVtkObject = chartViewNodes->GetNextItemAsObject();
  vtkMRMLChartViewNode* chartViewNode = vtkMRMLChartViewNode::SafeDownCast(chartViewNodeVtkObject);*/
  // Get the Chart View Node
  vtkMRMLChartViewNode* chartViewNode = vtkMRMLChartViewNode::SafeDownCast(mrmlScene->GetNextNodeByClass("vtkMRMLChartViewNode") );

  // Create an Array Node and add the eigen value
  VectorType eigenvalue = statModel->GetPCAVarianceVector();
  vtkMRMLDoubleArrayNode* doubleArrayNode = vtkMRMLDoubleArrayNode::SafeDownCast(mrmlScene->AddNode(vtkMRMLDoubleArrayNode::New()));
  vtkDoubleArray* a = doubleArrayNode->GetArray();
  int nbPrincipalComponent = statModel->GetNumberOfPrincipalComponents();
  a->SetNumberOfTuples(nbPrincipalComponent);
  for (int i = 0;i<nbPrincipalComponent;i++){
    std::cout<<"i "<<i<<std::endl;
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
  
  // Tell the Chart View which Chart to display
  chartViewNode->SetChartNodeID(chartNode->GetID());
}

void vtkSlicerLoadableSSMBuildingLogic::ObserveMRMLScene()
{
  this->Superclass::ObserveMRMLScene();
}


