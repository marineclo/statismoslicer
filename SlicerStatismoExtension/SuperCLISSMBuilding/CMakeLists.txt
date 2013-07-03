
#-----------------------------------------------------------------------------
# Third party libraries
find_package(Statismo REQUIRED)

link_directories(${STATISMO_LIBRARY_DIR})

#-----------------------------------------------------------------------------
set(MODULE_NAME SuperCLISSMBuilding)

#-----------------------------------------------------------------------------
set(MODULE_INCLUDE_DIRECTORIES
  ${STATISMO_INCLUDE_DIRS}
  )

# Variable STATISMO_DIR is defined in 'lib/cmake/statismo-0.8/statismo-config.cmake'
set(MODULE_SRCS
  ${STATISMO_DIR}/include/Representers/VTK/vtkPolyDataRepresenter.h
  )

set(MODULE_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  ${STATISMO_LIBRARIES}
  ${VTK_LIBRARIES} ${VTK_HDF5_LIBRARIES}
  )
message("STATISMO_LIBRARIES: ${STATISMO_LIBRARIES}")
message("VTK_LIBRARIES: ${VTK_LIBRARIES}")
#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${MODULE_NAME}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  ADDITIONAL_SRCS ${MODULE_SRCS}
  EXECUTABLE_ONLY
  )

#-----------------------------------------------------------------------------
#if(BUILD_TESTING)
#  add_subdirectory(Testing)
#endif()

