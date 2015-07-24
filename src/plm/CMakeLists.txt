add_executable(ad_plm ad_plm.c ad_lib_plm.c ad_lib_plm.h)
target_link_libraries(ad_plm ad)

add_dependencies(check ad_plm)

# Add a test for fms
add_test(NAME plm COMMAND utest.sh ${CMAKE_CURRENT_BINARY_DIR}
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
         
