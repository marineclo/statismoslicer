#ifndef __qSlicerLoadableSSMBuildingSelectParametersWidget_h
#define __qSlicerLoadableSSMBuildingSelectParametersWidget_h

// Qt includes
#include <QWidget>

// FooBar Widgets includes
#include "qSlicerLoadableSSMBuildingModuleWidgetsExport.h"

class qSlicerLoadableSSMBuildingSelectParametersWidgetPrivate;

/// \ingroup Slicer_QtModules_LoadableSSMBuilding
class Q_SLICER_MODULE_LOADABLESSMBUILDING_WIDGETS_EXPORT qSlicerLoadableSSMBuildingSelectParametersWidget
  : public QWidget
{
  Q_OBJECT
public:
  typedef QWidget Superclass;
  qSlicerLoadableSSMBuildingSelectParametersWidget(QWidget *parent=0);
  virtual ~qSlicerLoadableSSMBuildingSelectParametersWidget();

protected slots:
  void onSelect();
  void onSelectInputModel();

protected:
  QScopedPointer<qSlicerLoadableSSMBuildingSelectParametersWidgetPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerLoadableSSMBuildingSelectParametersWidget);
  Q_DISABLE_COPY(qSlicerLoadableSSMBuildingSelectParametersWidget);
};

#endif
