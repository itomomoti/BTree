# set(CMAKE_VERBOSE_MAKEFILE 1)
cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) ## on g++ this ensures: -std=c++11 and not -std=gnu++11



#### Build setting
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  #### Set a default build type if none was specified
  if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  endif()
  message("Building in ${CMAKE_BUILD_TYPE} mode")

  set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -p -Wconversion")
  set(CMAKE_CXX_FLAGS_RELEASE "-ggdb -Ofast -fstrict-aliasing -DNDEBUG -march=native")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")
endif()



include_directories(${CMAKE_SOURCE_DIR}/module/Basics)



#### set runtime directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})



#### Generate executable file
add_executable(PSumWithValue_Test ${CMAKE_CURRENT_SOURCE_DIR}/PSumWithValue_Test.cpp ${CMAKE})

# add_executable(BitsUtil_Benchmark ${CMAKE_CURRENT_SOURCE_DIR}/BitsUtil_Benchmark.cpp ${CMAKE})

#### Link libraries
target_link_libraries(PSumWithValue_Test Basics gtest gtest_main pthread)

# target_link_libraries(BitsUtil_Benchmark Basics)
