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

#ifndef __qSlicerSSMViewerFooBarWidget_h
#define __qSlicerSSMViewerFooBarWidget_h

// Qt includes
#include <QWidget>

// FooBar Widgets includes
#include "qSlicerSSMViewerModuleWidgetsExport.h"

class qSlicerSSMViewerFooBarWidgetPrivate;

/// \ingroup Slicer_QtModules_SSMViewer
class Q_SLICER_MODULE_SSMVIEWER_WIDGETS_EXPORT qSlicerSSMViewerFooBarWidget
  : public QWidget
{
  Q_OBJECT
public:
  typedef QWidget Superclass;
  qSlicerSSMViewerFooBarWidget(QWidget *parent=0);
  virtual ~qSlicerSSMViewerFooBarWidget();

protected slots:

protected:
  QScopedPointer<qSlicerSSMViewerFooBarWidgetPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerSSMViewerFooBarWidget);
  Q_DISABLE_COPY(qSlicerSSMViewerFooBarWidget);
};

#endif
