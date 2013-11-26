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

//Statismo includes
#include "itkMeshRepresenter.h"
#include "statismo_ITK/itkStatisticalModel.h"
#include "Representers/VTK/vtkPolyDataRepresenter.h"

// For VTK model
typedef vtkPolyDataRepresenter VtkRepresenterType;
typedef statismo::StatisticalModel<VtkRepresenterType> VtkStatisticalModelType;
using std::auto_ptr;

//Fot ITK model
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
  void displayEigenSpectrum(unsigned int nbPrincipalComponent);
  ItkStatisticalModelType::Pointer itkModel;
  //auto_ptr<VtkStatisticalModelType> vtkModel;
  VtkStatisticalModelType* vtkModel;
  std::vector<std::string> nameITK;
  int indexNode;

public slots:
  void onSelect();
  void onSelectInputModel();
  void applyModel();

protected:
  QScopedPointer<qSlicerDisplaySSMModuleWidgetPrivate> d_ptr;
  
  virtual void setup();
  virtual void setMRMLScene(vtkMRMLScene *scene);

private:
  Q_DECLARE_PRIVATE(qSlicerDisplaySSMModuleWidget);
  Q_DISABLE_COPY(qSlicerDisplaySSMModuleWidget);
};

#endif
