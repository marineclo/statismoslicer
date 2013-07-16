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

// STD includes
#include <cassert>

#include "vtkSlicerModelsLogic.h"
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
void vtkSlicerLoadableSSMBuildingLogic::DisplaySampleModel(vtkPolyData* polydata)
{
  /*vtkSlicerModelsLogic* modelsLogic = vtkSlicerModelsLogic::New();
  vtkMRMLModelNode* modelNode = modelsLogic->AddModel(polydata);
  if (modelNode == NULL){
    std::cout<<"empty"<<std::endl;
  }*/
  std::cout<<"output "<<this->GetMRMLScene()<<std::endl;
  vtkMRMLScene* mrmlScene = vtkMRMLScene::New();
  vtkMRMLModelNode* sampleNode = vtkMRMLModelNode::New();
  sampleNode->SetScene(mrmlScene);
  sampleNode->SetName("Sample");
  sampleNode->SetAndObservePolyData(polydata);
  
  vtkMRMLModelDisplayNode* modelDisplayNode = vtkMRMLModelDisplayNode::New();
  modelDisplayNode->SetColor(0,1,0); // green
  modelDisplayNode->SetScene(mrmlScene);
  mrmlScene->AddNode(modelDisplayNode);
  sampleNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());
  
  modelDisplayNode->SetInputPolyData(sampleNode->GetPolyData());
  mrmlScene->AddNode(sampleNode);
  
  std::cout<<"output1 "<<this->GetMRMLScene()<<std::endl;
}

