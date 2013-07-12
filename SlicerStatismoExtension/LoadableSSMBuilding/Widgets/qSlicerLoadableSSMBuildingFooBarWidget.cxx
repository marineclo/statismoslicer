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
#include "qSlicerLoadableSSMBuildingFooBarWidget.h"
#include "ui_qSlicerLoadableSSMBuildingFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_LoadableSSMBuilding
class qSlicerLoadableSSMBuildingFooBarWidgetPrivate
  : public Ui_qSlicerLoadableSSMBuildingFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerLoadableSSMBuildingFooBarWidget);
protected:
  qSlicerLoadableSSMBuildingFooBarWidget* const q_ptr;

public:
  qSlicerLoadableSSMBuildingFooBarWidgetPrivate(
    qSlicerLoadableSSMBuildingFooBarWidget& object);
  virtual void setupUi(qSlicerLoadableSSMBuildingFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerLoadableSSMBuildingFooBarWidgetPrivate
::qSlicerLoadableSSMBuildingFooBarWidgetPrivate(
  qSlicerLoadableSSMBuildingFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerLoadableSSMBuildingFooBarWidgetPrivate
::setupUi(qSlicerLoadableSSMBuildingFooBarWidget* widget)
{
  this->Ui_qSlicerLoadableSSMBuildingFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerLoadableSSMBuildingFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerLoadableSSMBuildingFooBarWidget
::qSlicerLoadableSSMBuildingFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerLoadableSSMBuildingFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerLoadableSSMBuildingFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerLoadableSSMBuildingFooBarWidget
::~qSlicerLoadableSSMBuildingFooBarWidget()
{
}
