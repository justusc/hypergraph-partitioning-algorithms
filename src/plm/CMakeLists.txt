add_library(libplm libplm.c)

add_executable(plm plm.c)
target_link_libraries(plm libplm libhgpa)

add_dependencies(check plm)

# Add a test for fms
add_test(NAME plm COMMAND utest.sh ${CMAKE_CURRENT_BINARY_DIR}
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
         
