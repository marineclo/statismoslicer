
#-----------------------------------------------------------------------------
# Third party libraries
find_package(Statismo REQUIRED)

link_directories(${STATISMO_LIBRARY_DIR})

#-----------------------------------------------------------------------------
project(qSlicer${MODULE_NAME}ModuleWidgets)

set(KIT ${PROJECT_NAME})

set(${KIT}_EXPORT_DIRECTIVE "Q_SLICER_MODULE_${MODULE_NAME_UPPER}_WIDGETS_EXPORT")

set(${KIT}_INCLUDE_DIRECTORIES
  ${STATISMO_INCLUDE_DIRS}
  )

set(${KIT}_SRCS
  qSlicer${MODULE_NAME}FooBarWidget.cxx
  qSlicer${MODULE_NAME}FooBarWidget.h
  qSlicer${MODULE_NAME}SelectParametersWidget.cxx
  qSlicer${MODULE_NAME}SelectParametersWidget.h
  )

set(${KIT}_MOC_SRCS
  qSlicer${MODULE_NAME}FooBarWidget.h
  qSlicer${MODULE_NAME}SelectParametersWidget.h
  )

set(${KIT}_UI_SRCS
  ../Resources/UI/qSlicer${MODULE_NAME}FooBarWidget.ui
  ../Resources/UI/qSlicer${MODULE_NAME}SelectParametersWidget.ui
  )

set(${KIT}_RESOURCES
  ../Resources/UI/qSlicer${MODULE_NAME}FooBarWidget.ui
  ../Resources/UI/qSlicer${MODULE_NAME}SelectParametersWidget.ui
  )

set(${KIT}_TARGET_LIBRARIES
  vtkSlicer${MODULE_NAME}ModuleLogic
  ${ITK_LIBRARIES}
  ${STATISMO_LIBRARIES}
  ${VTK_LIBRARIES} ${VTK_HDF5_LIBRARIES}
  )

#-----------------------------------------------------------------------------
SlicerMacroBuildModuleWidgets(
  NAME ${KIT}
  EXPORT_DIRECTIVE ${${KIT}_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${${KIT}_INCLUDE_DIRECTORIES}
  SRCS ${${KIT}_SRCS}
  MOC_SRCS ${${KIT}_MOC_SRCS}
  UI_SRCS ${${KIT}_UI_SRCS}
  TARGET_LIBRARIES ${${KIT}_TARGET_LIBRARIES}
  RESOURCES ${${KIT}_RESOURCES}
  WRAP_PYTHONQT
  )

