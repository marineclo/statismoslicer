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

// .NAME vtkSlicerSSMViewerLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerSSMViewerLogic_h
#define __vtkSlicerSSMViewerLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"
class vtkSlicerModelsLogic;

// VTK includes
class vtkPolyData;

#include "vtkSlicerSSMViewerModuleLogicExport.h"

/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_SSMVIEWER_MODULE_LOGIC_EXPORT vtkSlicerSSMViewerLogic :
  public vtkSlicerModuleLogic
{
public:

  static vtkSlicerSSMViewerLogic *New();
  vtkTypeMacro(vtkSlicerSSMViewerLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);
  void DisplaySampleModel(vtkPolyData* polydata, vtkMRMLScene * mrmlScene);

protected:
  vtkSlicerSSMViewerLogic();
  virtual ~vtkSlicerSSMViewerLogic();

  virtual void ObserveMRMLScene();
  virtual void SetMRMLSceneInternal(vtkMRMLScene* newScene);
  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  virtual void RegisterNodes();
  virtual void UpdateFromMRMLScene();
  virtual void OnMRMLSceneNodeAdded(vtkMRMLNode* node);
  virtual void OnMRMLSceneNodeRemoved(vtkMRMLNode* node);
private:

  vtkSlicerSSMViewerLogic(const vtkSlicerSSMViewerLogic&); // Not implemented
  void operator=(const vtkSlicerSSMViewerLogic&);               // Not implemented
};

#endif
