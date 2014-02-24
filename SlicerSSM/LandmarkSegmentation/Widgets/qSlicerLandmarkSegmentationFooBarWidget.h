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

#ifndef __qSlicerLandmarkSegmentationFooBarWidget_h
#define __qSlicerLandmarkSegmentationFooBarWidget_h

// Qt includes
#include <QWidget>

// FooBar Widgets includes
#include "qSlicerLandmarkSegmentationModuleWidgetsExport.h"

class qSlicerLandmarkSegmentationFooBarWidgetPrivate;

/// \ingroup Slicer_QtModules_LandmarkSegmentation
class Q_SLICER_MODULE_LANDMARKSEGMENTATION_WIDGETS_EXPORT qSlicerLandmarkSegmentationFooBarWidget
  : public QWidget
{
  Q_OBJECT
public:
  typedef QWidget Superclass;
  qSlicerLandmarkSegmentationFooBarWidget(QWidget *parent=0);
  virtual ~qSlicerLandmarkSegmentationFooBarWidget();

protected slots:

protected:
  QScopedPointer<qSlicerLandmarkSegmentationFooBarWidgetPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerLandmarkSegmentationFooBarWidget);
  Q_DISABLE_COPY(qSlicerLandmarkSegmentationFooBarWidget);
};

#endif
