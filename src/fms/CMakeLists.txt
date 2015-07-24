add_executable(ad_fms ad_fms.c ad_lib_fms.c ad_lib_fms.h)
target_link_libraries(ad_fms ad)


add_dependencies(check ad_fms)

# Add a test for fms
add_test(NAME fms COMMAND utest.sh ${CMAKE_CURRENT_BINARY_DIR}
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
         
