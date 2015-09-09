add_library(libfms libfms.c)

add_executable(fms fms.c)
target_link_libraries(fms libfms libhgpa)


add_dependencies(check fms)

# Add a test for fms
add_test(NAME fms COMMAND utest.sh ${CMAKE_CURRENT_BINARY_DIR}
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
         
