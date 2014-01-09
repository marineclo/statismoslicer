
# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

# Sanity checks
if(DEFINED Statismo_DIR AND NOT EXISTS ${Statismo_DIR})
  message(FATAL_ERROR "Statismo_DIR variable is defined but corresponds to non-existing directory")
endif()

# Set dependency list
set(Statismo_DEPENDENCIES "")

# Include dependent projects if any
SlicerMacroCheckExternalProjectDependency(Statismo)
set(proj Statismo)

if(NOT DEFINED ${proj}_DIR)
  #message(STATUS "${__indent}Adding project ${proj}")

  set(EXTERNAL_PROJECT_OPTIONAL_ARGS)

  # Set CMake OSX variable to pass down the external project
  if(APPLE)
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
      -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
      -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
      -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET})
  endif()

  option(CHOOSE_GIT_PROTOCOL "Choose to use git or https protocol" ON)
  if(CHOOSE_GIT_PROTOCOL)
    set(git_protocol "git")
  else(CHOOSE_GIT_PROTOCOL)
    set(git_protocol "https")
  endif()
  
  #if(NOT DEFINED git_protocol)
  #  set(git_protocol "git")
  #endif()

  ExternalProject_Add(${proj}
    GIT_REPOSITORY "${git_protocol}://github.com/statismo/statismo.git"
    #GIT_REPOSITORY "${git_protocol}://github.com/arnaudgelas/statismo.git"
    GIT_TAG "b9ad869b4df1f258f59c4cc31ef2f76c69b965ec"
    #GIT_TAG "BuildSystem"
    SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj}
    BINARY_DIR ${proj}-build
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      #-DCMAKE_CXX_FLAGS:STRING=${ep_common_cxx_flags}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      #-DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/${proj}-install
      #-DBUILD_TESTING:BOOL=OFF
      ${EXTERNAL_PROJECT_OPTIONAL_ARGS}
    DEPENDS
      ${${proj}_DEPENDENCIES}
    )
    
  #set(${proj}_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
  #set(${proj}_DIR ${CMAKE_BINARY_DIR}/${proj}-install)
  set(${proj}_DIR ${CMAKE_BINARY_DIR}/${proj}-install/lib/cmake/statismo-0.81)
  
else()
  # The project is provided using <proj>_DIR, nevertheless since other project may depend on <proj>,
  # let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${${proj}_DEPENDENCIES}")
endif()

