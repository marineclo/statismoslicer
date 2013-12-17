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

#ifndef __qSlicerLandmarkSegmentationModuleWidget_h
#define __qSlicerLandmarkSegmentationModuleWidget_h

// SlicerQt includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerLandmarkSegmentationModuleExport.h"

class qSlicerLandmarkSegmentationModuleWidgetPrivate;
class vtkMRMLNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_LANDMARKSEGMENTATION_EXPORT qSlicerLandmarkSegmentationModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerLandmarkSegmentationModuleWidget(QWidget *parent=0);
  virtual ~qSlicerLandmarkSegmentationModuleWidget();

public slots:
  void setModel();
  //void setCTscan();
  void apply();


protected:
  QScopedPointer<qSlicerLandmarkSegmentationModuleWidgetPrivate> d_ptr;
  
  virtual void setup();
  virtual void setMRMLScene(vtkMRMLScene *scene);

private:
  Q_DECLARE_PRIVATE(qSlicerLandmarkSegmentationModuleWidget);
  Q_DISABLE_COPY(qSlicerLandmarkSegmentationModuleWidget);
};

#endif
