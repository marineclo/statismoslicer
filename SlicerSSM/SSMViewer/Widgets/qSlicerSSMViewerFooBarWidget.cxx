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
#include "qSlicerSSMViewerFooBarWidget.h"
#include "ui_qSlicerSSMViewerFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_SSMViewer
class qSlicerSSMViewerFooBarWidgetPrivate
  : public Ui_qSlicerSSMViewerFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerSSMViewerFooBarWidget);
protected:
  qSlicerSSMViewerFooBarWidget* const q_ptr;

public:
  qSlicerSSMViewerFooBarWidgetPrivate(
    qSlicerSSMViewerFooBarWidget& object);
  virtual void setupUi(qSlicerSSMViewerFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerSSMViewerFooBarWidgetPrivate
::qSlicerSSMViewerFooBarWidgetPrivate(
  qSlicerSSMViewerFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerSSMViewerFooBarWidgetPrivate
::setupUi(qSlicerSSMViewerFooBarWidget* widget)
{
  this->Ui_qSlicerSSMViewerFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerSSMViewerFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerSSMViewerFooBarWidget
::qSlicerSSMViewerFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerSSMViewerFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerSSMViewerFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerSSMViewerFooBarWidget
::~qSlicerSSMViewerFooBarWidget()
{
}
