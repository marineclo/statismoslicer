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

#ifndef __qSlicerLoadableSSMBuildingModuleWidget_h
#define __qSlicerLoadableSSMBuildingModuleWidget_h

// SlicerQt includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerLoadableSSMBuildingModuleExport.h"
#include "itkPoint.h"

typedef std::vector<std::string> StringVectorType;
typedef itk::Point<double, 3> PointType;

class qSlicerLoadableSSMBuildingModuleWidgetPrivate;
class vtkMRMLNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_LOADABLESSMBUILDING_EXPORT qSlicerLoadableSSMBuildingModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerLoadableSSMBuildingModuleWidget(QWidget *parent=0);
  virtual ~qSlicerLoadableSSMBuildingModuleWidget();
  void displayEigenSpectrum();
  std::vector<PointType> readLandmarks(const std::string& filename);
  int getdir (std::string dir, StringVectorType &files, std::string nameImage, const std::string& extension, std::string &ref, std::string numRef);
  int getdir2 (std::string dir, StringVectorType &files, std::string nameImage, const std::string& extension);

public slots:
  void onSelect();
  void onSelectInputModel();
  void onPerformedRegistration();
  void onBuildPartiallyFixedModel();
  void onShapeModelFitting();

protected:
  QScopedPointer<qSlicerLoadableSSMBuildingModuleWidgetPrivate> d_ptr;
  
  virtual void setup();
  virtual void setMRMLScene(vtkMRMLScene *scene);

private:
  Q_DECLARE_PRIVATE(qSlicerLoadableSSMBuildingModuleWidget);
  Q_DISABLE_COPY(qSlicerLoadableSSMBuildingModuleWidget);
};

#endif
