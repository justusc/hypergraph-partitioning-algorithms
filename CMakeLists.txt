cmake_minimum_required(VERSION 3.1.0)
project(HGPA C)

include(CMakeDependentOption)

get_property(SUPPORTS_SHARED GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS)
cmake_dependent_option(BUILD_SHARED_LIBS "Enable shared libraries" ON
      "SUPPORTS_SHARED;ENABLE_SHARED_LIBRARIES" OFF)
option(ENABLE_TRACE "Enable partition trace" OFF)

if(ENABLE_TRACE)
  add_definitions(-DAD_TRACE=1)
endif()

# Add c99 Standard flag
if (CMAKE_VERSION VERSION_LESS "3.1")
  set(CMAKE_C_FLAGS "--std=c99 ${CMAKE_C_FLAGS}")
else ()
  set(CMAKE_C_STANDARD 99)
  set(CMAKE_C_EXTENSIONS OFF)
endif ()

include_directories(src/share)

# Add check target
enable_testing()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} -V)

include_directories(include)

add_subdirectory(src/share)
add_subdirectory(src/fms)
add_subdirectory(src/pfm)
add_subdirectory(src/plm)
add_subdirectory(src/sa)