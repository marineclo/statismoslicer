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

#ifndef __qSlicerDisplaySSMModule_h
#define __qSlicerDisplaySSMModule_h

// SlicerQt includes
#include "qSlicerLoadableModule.h"

#include "qSlicerDisplaySSMModuleExport.h"

class qSlicerDisplaySSMModulePrivate;

/// \ingroup Slicer_QtModules_ExtensionTemplate
class Q_SLICER_QTMODULES_DISPLAYSSM_EXPORT
qSlicerDisplaySSMModule
  : public qSlicerLoadableModule
{
  Q_OBJECT
  Q_INTERFACES(qSlicerLoadableModule);

public:

  typedef qSlicerLoadableModule Superclass;
  explicit qSlicerDisplaySSMModule(QObject *parent=0);
  virtual ~qSlicerDisplaySSMModule();

  qSlicerGetTitleMacro(QTMODULE_TITLE);

  virtual QString helpText()const;
  virtual QString acknowledgementText()const;
  virtual QStringList contributors()const;

  virtual QIcon icon()const;

  virtual QStringList categories()const;
  virtual QStringList dependencies() const;

//public slots:
  virtual void setMRMLScene(vtkMRMLScene* mrmlScene);

protected:

  /// Initialize the module. Register the volumes reader/writer
  virtual void setup();

  /// Create and return the widget representation associated to this module
  virtual qSlicerAbstractModuleRepresentation * createWidgetRepresentation();

  /// Create and return the logic associated to this module
  virtual vtkMRMLAbstractLogic* createLogic();

protected:
  QScopedPointer<qSlicerDisplaySSMModulePrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerDisplaySSMModule);
  Q_DISABLE_COPY(qSlicerDisplaySSMModule);

};

#endif
