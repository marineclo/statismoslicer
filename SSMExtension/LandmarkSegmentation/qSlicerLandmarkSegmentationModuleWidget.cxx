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
#include "qSlicerLandmarkSegmentationModuleWidget.h"
#include "ui_qSlicerLandmarkSegmentationModuleWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerLandmarkSegmentationModuleWidgetPrivate: public Ui_qSlicerLandmarkSegmentationModuleWidget
{
public:
  qSlicerLandmarkSegmentationModuleWidgetPrivate();
};

//-----------------------------------------------------------------------------
// qSlicerLandmarkSegmentationModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerLandmarkSegmentationModuleWidgetPrivate::qSlicerLandmarkSegmentationModuleWidgetPrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerLandmarkSegmentationModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerLandmarkSegmentationModuleWidget::qSlicerLandmarkSegmentationModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerLandmarkSegmentationModuleWidgetPrivate )
{
}

//-----------------------------------------------------------------------------
qSlicerLandmarkSegmentationModuleWidget::~qSlicerLandmarkSegmentationModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerLandmarkSegmentationModuleWidget::setup()
{
  Q_D(qSlicerLandmarkSegmentationModuleWidget);
  d->setupUi(this);
  //this->Superclass::setup();
}

//-----------------------------------------------------------------------------
void qSlicerLandmarkSegmentationModuleWidget::setMRMLScene(vtkMRMLScene* mrmlScene)
{
  this->qSlicerAbstractModuleWidget::setMRMLScene(mrmlScene);
}


void qSlicerLandmarkSegmentationModuleWidget::setRefFid(){

  Q_D(qSlicerLandmarkSegmentationModuleWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select reference fiducial", QString());
  d->fixedLandmarksName->setText(inputFile);
}

void qSlicerLandmarkSegmentationModuleWidget::setModel(){
  Q_D(qSlicerLandmarkSegmentationModuleWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select model", QString());
  d->modelName->setText(inputFile);
}

void qSlicerLandmarkSegmentationModuleWidget::setInpuFid(){
  Q_D(qSlicerLandmarkSegmentationModuleWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select input Fiducial", QString());
  d->movingLandmarksName->setText(inputFile);
}

void qSlicerLandmarkSegmentationModuleWidget::setCTscan(){
  Q_D(qSlicerLandmarkSegmentationModuleWidget);
  QString inputFile = QFileDialog::getOpenFileName(this, "Select CT scan", QString());
  d->targetName->setText(inputFile);
}

void qSlicerLandmarkSegmentationModuleWidget::apply(){

}
