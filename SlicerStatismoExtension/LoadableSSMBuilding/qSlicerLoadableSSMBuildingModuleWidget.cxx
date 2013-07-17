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

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerLoadableSSMBuildingModuleWidgetPrivate: public Ui_qSlicerLoadableSSMBuildingModuleWidget
{
public:
  qSlicerLoadableSSMBuildingModuleWidgetPrivate();
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

}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::onSelect()
{
  /**** VTK Model ****/
  using std::auto_ptr;
  using namespace statismo;
  typedef vtkPolyDataRepresenter RepresenterType;
  typedef StatisticalModel<RepresenterType> StatisticalModelType;
  Q_D(qSlicerLoadableSSMBuildingModuleWidget);
  // Get the model name selected by the user
  QString modelQString = d->modelNamePath->text();
  std::string modelString = modelQString.toStdString(); 
  // Load the model
  auto_ptr<StatisticalModelType> statModel(StatisticalModelType::Load(modelString.c_str()));
  VectorType coefficients = VectorType::Zero(statModel->GetNumberOfPrincipalComponents());
  int pc = static_cast<int>(d->pcSlider->value());
  coefficients(pc) = d->stdSlider->value();
  vtkPolyData* samplePC = statModel->DrawSample(coefficients);
  
  std::cout<<"output "<<this->mrmlScene()<<std::endl;
  // Add polydata sample to the scene
  //vtkSlicerLoadableSSMBuildingLogic* moduleLogic = vtkSlicerLoadableSSMBuildingLogic::New();
  //moduleLogic->DisplaySampleModel(samplePC);
  vtkMRMLScene* mrmlScene = this->mrmlScene();
  vtkMRMLModelNode* sampleNode = vtkMRMLModelNode::New();
  sampleNode->SetScene(mrmlScene);
  sampleNode->SetName("Sample");
  sampleNode->SetAndObservePolyData(samplePC);
  
  vtkMRMLModelDisplayNode* modelDisplayNode = vtkMRMLModelDisplayNode::New();
  modelDisplayNode->SetColor(0,1,0); // green
  modelDisplayNode->SetScene(mrmlScene);
  mrmlScene->AddNode(modelDisplayNode);
  sampleNode->SetAndObserveDisplayNodeID(modelDisplayNode->GetID());
  
  modelDisplayNode->SetInputPolyData(sampleNode->GetPolyData());
  mrmlScene->AddNode(sampleNode);
  
}

//-----------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingModuleWidget::setMRMLScene(vtkMRMLScene* mrmlScene)
{
  this->qSlicerAbstractModuleWidget::setMRMLScene(mrmlScene);
}

