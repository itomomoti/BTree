# set(CMAKE_VERBOSE_MAKEFILE 1)
cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) ## on g++ this ensures: -std=c++11 and not -std=gnu++11



#### Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()
message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -p -Wconversion")
set(CMAKE_CXX_FLAGS_RELEASE "-ggdb -Ofast -fstrict-aliasing -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

project(BTree)



#### Modules
add_subdirectory(${PROJECT_SOURCE_DIR}/module/Basics)
include_directories(${PROJECT_SOURCE_DIR}/module/Basics)

include_directories(${CMAKE_SOURCE_DIR})

#add_library(${PROJECT_NAME} STATIC
#  hoge.cpp
#)



#### TEST (to enable this add "-DTESTING=1" when running cmake)
if(TESTING)
  add_subdirectory(test)
  enable_testing()
endif(TESTING)



#### Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  set(DOXYGEN_CONF_FILE "doxygen.conf")
  set(DOXYGEN_CONF_PATH ${CMAKE_CURRENT_BINARY_DIR}/${DOXYGEN_CONF_FILE})
  set(DOXYGEN_SOURCE_DIR ${CMAKE_SOURCE_DIR})
  set(DOXYGEN_TARGET "doc")

  configure_file(${DOXYGEN_CONF_FILE}.in ${DOXYGEN_CONF_PATH})

  add_custom_target(${DOXYGEN_TARGET}
                    ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONF_PATH}
		    DEPENDS ${DOXYGEN_CONF_PATH})
  message("\"$ make doc\" to run Doxygen")
else(DOXYGEN_FOUND)
  message(WARNING "doxygen binary couldn't be found")
endif(DOXYGEN_FOUND)
