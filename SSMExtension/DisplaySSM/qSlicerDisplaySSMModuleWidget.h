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

#ifndef __qSlicerDisplaySSMModuleWidget_h
#define __qSlicerDisplaySSMModuleWidget_h

// SlicerQt includes
#include "qSlicerAbstractModuleWidget.h"

#include "qSlicerDisplaySSMModuleExport.h"

#include "itkMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"

const unsigned Dimensions = 3;
typedef itk::MeshRepresenter<float, Dimensions> ItkRepresenterType;
typedef itk::StatisticalModel<ItkRepresenterType> ItkStatisticalModelType;
typedef vnl_vector<statismo::ScalarType> itkVectorType;


class qSlicerDisplaySSMModuleWidgetPrivate;
class vtkMRMLNode;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_DISPLAYSSM_EXPORT qSlicerDisplaySSMModuleWidget :
  public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:

  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerDisplaySSMModuleWidget(QWidget *parent=0);
  virtual ~qSlicerDisplaySSMModuleWidget();
  //void displayEigenSpectrum(ItkStatisticalModelType* statModel);
  void displayEigenSpectrum(itkVectorType itkEigenvalue, unsigned int nbPrincipalComponent);
  ItkStatisticalModelType::Pointer modelITK;


public slots:
  void onSelect();
  void onSelectInputModel();

protected:
  QScopedPointer<qSlicerDisplaySSMModuleWidgetPrivate> d_ptr;
  
  virtual void setup();
  virtual void setMRMLScene(vtkMRMLScene *scene);

private:
  Q_DECLARE_PRIVATE(qSlicerDisplaySSMModuleWidget);
  Q_DISABLE_COPY(qSlicerDisplaySSMModuleWidget);
};

#endif