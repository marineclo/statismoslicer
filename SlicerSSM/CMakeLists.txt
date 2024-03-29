cmake_minimum_required(VERSION 2.8.9)

#-----------------------------------------------------------------------------
project(SlicerSSM)

#set(EXTENSION_NAME SlicerSSM)

#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/4.3/Extensions/SlicerSSM")
set(EXTENSION_CATEGORY "Shape Analysis")
set(EXTENSION_CONTRIBUTORS "Marine Clogenson (EPFL), Jean-Christophe Fillon-Robin (Kitware)")
set(EXTENSION_DESCRIPTION "The SlicerSSM extension allow working with Statistical Shape Model (SSM). This extension uses Statismo, a framework for Statistical Shape Model building.")
set(EXTENSION_ICONURL "/Users/Marine/EPFL/Documents/Rapport/paper/screenshotSlicerView.png")
#set(EXTENSION_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/SuperBuildLoadableExtensionTemplate/Resources/Icons/SuperBuildLoadableExtensionTemplate.png?revision=19437&view=co")
#set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/b/b3/Slicer-r19441-SuperBuildLoadableExtensionTemplate-screenshot.png")
set(EXTENSION_DEPENDS "HDF5 Statismo") # Specified as a space separated list or 'NA' if any
set(EXTENSION_BUILD_SUBDIRECTORY inner-build)

set(SUPERBUILD_TOPLEVEL_PROJECT inner)

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})
mark_as_superbuild(Slicer_DIR)

find_package(Git REQUIRED)
mark_as_superbuild(GIT_EXECUTABLE)

#-----------------------------------------------------------------------------
option(${EXTENSION_NAME}_SUPERBUILD "Build ${EXTENSION_NAME} and the projects it depends on via SuperBuild.cmake." ON)
mark_as_advanced(${EXTENSION_NAME}_SUPERBUILD)
if(${EXTENSION_NAME}_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
endif()

#-----------------------------------------------------------------------------
set(Slicer_EXTENSION_EXTERNAL_PROJECT_LIBRARY_DIRS ${HDF5_DIR}/lib)

#-----------------------------------------------------------------------------
add_subdirectory(SSMViewer)
#add_subdirectory(LandmarkSegmentation)

#-----------------------------------------------------------------------------
set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/")
if(NOT APPLE)
  install(FILES ${HDF5_LIBRARIES} DESTINATION ${Slicer_INSTALL_ROOT}${Slicer_BUNDLE_EXTENSIONS_LOCATION}${Slicer_LIB_DIR} COMPONENT RuntimeLibraries)
endif()
include(${Slicer_EXTENSION_CPACK})
