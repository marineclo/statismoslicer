cmake_minimum_required(VERSION 2.8.9)

#find_package(Slicer COMPONENTS ConfigurePrerequisites REQUIRED)

project(SlicerStatismoExtension)

#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/SlicerStatismoExtension")
set(EXTENSION_CATEGORY "Examples")
set(EXTENSION_CONTRIBUTORS "Marine Clogenson (EPFL)")
set(EXTENSION_DESCRIPTION "This is an extension bundling a CLI module, to generate an atlas, built using SuperBuild to build a dependent project named Statismo")
#set(EXTENSION_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/SuperBuildLoadableExtensionTemplate/Resources/Icons/SuperBuildLoadableExtensionTemplate.png?revision=19437&view=co")
#set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/b/b3/Slicer-r19441-SuperBuildLoadableExtensionTemplate-screenshot.png")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated list or 'NA' if any
set(EXTENSION_BUILD_SUBDIRECTORY inner-build)

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

find_package(Git REQUIRED)

#-----------------------------------------------------------------------------
option(${EXTENSION_NAME}_SUPERBUILD "Build ${EXTENSION_NAME} and the projects it depends on via SuperBuild.cmake." ON)
mark_as_advanced(${EXTENSION_NAME}_SUPERBUILD)
if(${EXTENSION_NAME}_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
endif()

#-----------------------------------------------------------------------------
add_subdirectory(SuperCLISSMBuilding)
add_subdirectory(ScriptedSSMBuilding)

#-----------------------------------------------------------------------------
set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/")
#set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${Foo_DIR};Foo;RuntimeLibraries;/")
INCLUDE(${Slicer_EXTENSION_CPACK})
