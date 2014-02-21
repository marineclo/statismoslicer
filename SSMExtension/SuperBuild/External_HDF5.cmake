
set(proj HDF5)

# Set dependency list
set(${proj}_DEPENDS "")

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  message(FATAL_ERROR "Enabling ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj} is not supported !")
endif()

# Sanity checks
if(DEFINED Foo_DIR AND NOT EXISTS ${Foo_DIR})
  message(FATAL_ERROR "Foo_DIR variable is defined but corresponds to non-existing directory")
endif()

if(NOT DEFINED ${proj}_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})

  if(NOT DEFINED git_protocol)
    set(git_protocol "git")
  endif()
  
  # it seems that on linux only the shared libraries work, while on windows static seems to work fine
  if(WIN32)
	  set(HDF5_BUILD_SHARED OFF)
	  set(HDF5_LIBRARIES debug hdf5d debug hdf5_cppd optimized hdf5; optimized hdf5_cpp)
	  set(HDF5_BUILD_TYPE "RELEASE" CACHE STRING "Build type (RELEASE DEBUG) of the hdf5 library" )
  else(WIN32)
	  set(HDF5_BUILD_SHARED ON)
	  set(HDF5_LIBRARIES debug hdf5 debug hdf5_cpp optimized hdf5; optimized hdf5_cpp)
	  set(HDF5_BUILD_TYPE "RELEASE")
  endif(WIN32)

  if(APPLE)
    set( HDF5_VERSION "1.8.10" )
  else(APPLE)
      set( HDF5_VERSION "1.8.10" )
  endif(APPLE)
  
  set(${proj}_DIR ${CMAKE_CURRENT_BINARY_DIR}/3rdParty/HDF5)
  
  if(UNIX)
    # Fix build error for the form 'H5detect.c:146:1: error: unknown type name ‘sigjmp_buf'
    # known to happen on 'gcc (Ubuntu/Linaro 4.8.1-10ubuntu9) 4.8.1'
    # Suggested fix from http://lists.boost.org/boost-build/2004/01/5512.php
    set(${proj}_CMAKE_C_FLAGS_INIT "-D_POSIX_SOURCE")
  endif()

  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    URL "ftp://www.hdfgroup.org/HDF5/releases/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz"
    #URL_MD5
    SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj}
    BINARY_DIR ${proj}-build
    CMAKE_CACHE_ARGS
      -DCMAKE_C_FLAGS_INIT:STRING=${${proj}_CMAKE_C_FLAGS_INIT}
      -DCMAKE_BUILD_TYPE:STRING=${HDF5_BUILD_TYPE}
	    -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=OFF
	    -DHDF5_BUILD_CPP_LIB:BOOL=ON
	    -DBUILD_SHARED_LIBS:BOOL=${HDF5_BUILD_SHARED}
	    -DHDF5_BUILD_TOOLS:BOOL=OFF
	    -DCMAKE_INSTALL_PREFIX:PATH=${${proj}_DIR}
	  INSTALL_DIR ${${proj}_DIR}
    DEPENDS
      ${${proj}_DEPENDS}
    )

  set(HDF5_INCLUDE_DIR ${${proj}_DIR}/include)
  set(HDF5_INCLUDE_DIR_CPP ${HDF5_INCLUDE_DIR}/cpp)
  set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR} ${HDF5_INCLUDE_DIR_CPP})

	set(HDF5_LIBRARY_DIR ${${proj}_DIR}/lib)
	ExternalProject_Message(${proj} "HDF5_LIBRARY_DIR:${HDF5_LIBRARY_DIR}")
	
	if(WIN32)
    set(HDF5_LIBRARY ${HDF5_LIBRARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}hdf5.lib)
    set(HDF5_CPP_LIBRARY ${HDF5_LIBRARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}hdf5_cpp.lib)
  else()
    set(HDF5_LIBRARY ${HDF5_LIBRARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}hdf5.so.1.8.10)
    set(HDF5_CPP_LIBRARY ${HDF5_LIBRARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}hdf5_cpp.so.1.8.10)
    if(APPLE)
      set(HDF5_LIBRARY ${HDF5_LIBRARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}hdf5.1.8.10.dylib)
      set(HDF5_CPP_LIBRARY ${HDF5_LIBRARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}hdf5_cpp.1.8.10.dylib)
    endif()
  endif()
  set(HDF5_LIBRARIES ${HDF5_LIBRARY} ${HDF5_CPP_LIBRARY})
  
  if(APPLE)
    ExternalProject_Add_Step(${proj} fix_rpath_hdf5
      COMMAND install_name_tool -id ${HDF5_LIBRARY} ${HDF5_LIBRARY}
      DEPENDEES install
      )
    ExternalProject_Add_Step(${proj} fix_rpath_hdf5_cpp
      COMMAND install_name_tool -id ${HDF5_CPP_LIBRARY} ${HDF5_CPP_LIBRARY}
      DEPENDEES install
      )
    # Resolve the absolute path of libhdf5.7.4.0.dylib in cd 3rdParty/HDF5/lib/   otool -L libhdf5_cpp.1.8.10.dylib
    # libhdf5_cpp.1.8.10.dylib link to libhdf5.1.8.10.dylib
    ExternalProject_Add_Step(${proj} link_hdf5.7.4.0_to_hdf5_cpp.1.8.10
      COMMAND install_name_tool -change libhdf5.7.4.0.dylib ${HDF5_LIBRARY} ${HDF5_CPP_LIBRARY}
      DEPENDEES install
      )
  endif()

else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDS})
endif()

mark_as_superbuild(${proj}_DIR:PATH)

mark_as_superbuild(VARS HDF5_INCLUDE_DIRS:STRING HDF5_LIBRARIES:STRING)

ExternalProject_Message(${proj} "HDF5_INCLUDE_DIR:${HDF5_INCLUDE_DIR}")
ExternalProject_Message(${proj} "HDF5_INCLUDE_DIR_CPP:${HDF5_INCLUDE_DIR_CPP}")
ExternalProject_Message(${proj} "HDF5_LIBRARY:${HDF5_LIBRARY}")
ExternalProject_Message(${proj} "HDF5_CPP_LIBRARY:${HDF5_CPP_LIBRARY}")

