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

// DisplaySSM Logic includes
#include "vtkSlicerDisplaySSMLogic.h"

// MRML includes

// VTK includes
#include <vtkNew.h>
#include "vtkMRMLModelNode.h"
#include "vtkMRMLModelDisplayNode.h"

// STD includes
#include <cassert>

#include "vtkSlicerModelsLogic.h"

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerDisplaySSMLogic);

//----------------------------------------------------------------------------
vtkSlicerDisplaySSMLogic::vtkSlicerDisplaySSMLogic()
{
}

//----------------------------------------------------------------------------
vtkSlicerDisplaySSMLogic::~vtkSlicerDisplaySSMLogic()
{
}

//----------------------------------------------------------------------------
void vtkSlicerDisplaySSMLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerDisplaySSMLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerDisplaySSMLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerDisplaySSMLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerDisplaySSMLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerDisplaySSMLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerDisplaySSMLogic::DisplaySampleModel(vtkPolyData* polydata, vtkMRMLScene * mrmlScene)
{
  this->SetMRMLScene(mrmlScene);
  vtkSlicerModelsLogic* modelsLogic = vtkSlicerModelsLogic::New();
  modelsLogic->SetMRMLScene(mrmlScene);
  vtkMRMLModelNode* modelNode = modelsLogic->AddModel(polydata);
  //modelNode->GetDisplayNode()->SetColor(1, 0, 0); 
  //modelNode->GetDisplayNode()->SetName();
  if (modelNode == NULL){
    std::cout<<"empty"<<std::endl;
  }
}

void vtkSlicerDisplaySSMLogic::ObserveMRMLScene()
{
  this->Superclass::ObserveMRMLScene();
}
