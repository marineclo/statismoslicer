/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

// FooBar Widgets includes
#include "qSlicerLandmarkSegmentationFooBarWidget.h"
#include "ui_qSlicerLandmarkSegmentationFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_LandmarkSegmentation
class qSlicerLandmarkSegmentationFooBarWidgetPrivate
  : public Ui_qSlicerLandmarkSegmentationFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerLandmarkSegmentationFooBarWidget);
protected:
  qSlicerLandmarkSegmentationFooBarWidget* const q_ptr;

public:
  qSlicerLandmarkSegmentationFooBarWidgetPrivate(
    qSlicerLandmarkSegmentationFooBarWidget& object);
  virtual void setupUi(qSlicerLandmarkSegmentationFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerLandmarkSegmentationFooBarWidgetPrivate
::qSlicerLandmarkSegmentationFooBarWidgetPrivate(
  qSlicerLandmarkSegmentationFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerLandmarkSegmentationFooBarWidgetPrivate
::setupUi(qSlicerLandmarkSegmentationFooBarWidget* widget)
{
  this->Ui_qSlicerLandmarkSegmentationFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerLandmarkSegmentationFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerLandmarkSegmentationFooBarWidget
::qSlicerLandmarkSegmentationFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerLandmarkSegmentationFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerLandmarkSegmentationFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerLandmarkSegmentationFooBarWidget
::~qSlicerLandmarkSegmentationFooBarWidget()
{
}
