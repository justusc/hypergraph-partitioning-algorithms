cmake_minimum_required(VERSION 3.1.0)
project(HGPA C)

# Add c99 Standard flag
if (CMAKE_VERSION VERSION_LESS "3.1")
  set(CMAKE_C_FLAGS "--std=c99 ${CMAKE_C_FLAGS}")
else ()
  set(CMAKE_C_STANDARD 99)
  set(CMAKE_C_EXTENSIONS OFF)
endif ()

include_directories(src/share)

add_subdirectory(src/share)
add_subdirectory(src/fms)
add_subdirectory(src/pfm)
add_subdirectory(src/plm)
add_subdirectory(src/sa)