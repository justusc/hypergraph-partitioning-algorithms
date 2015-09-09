add_library(libfms ad_lib_fms.c)

add_executable(ad_fms ad_fms.c)
target_link_libraries(ad_fms libfms libhgpa)


add_dependencies(check ad_fms)

# Add a test for fms
add_test(NAME fms COMMAND utest.sh ${CMAKE_CURRENT_BINARY_DIR}
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
         
