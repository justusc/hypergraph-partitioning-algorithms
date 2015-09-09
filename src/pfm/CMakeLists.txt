add_library(libpfm libpfm.c)

add_executable(pfm pfm.c)
target_link_libraries(pfm libpfm libhgpa)

add_dependencies(check pfm)

# Add a test for fms
add_test(NAME pfm COMMAND utest.sh ${CMAKE_CURRENT_BINARY_DIR}
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
         
