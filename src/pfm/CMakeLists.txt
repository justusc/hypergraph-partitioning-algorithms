add_library(libpfm ad_lib_pfm.c)

add_executable(ad_pfm ad_pfm.c)
target_link_libraries(ad_pfm libpfm libhgpa)

add_dependencies(check ad_pfm)

# Add a test for fms
add_test(NAME pfm COMMAND utest.sh ${CMAKE_CURRENT_BINARY_DIR}
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
         
