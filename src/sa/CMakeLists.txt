add_library(libsa libsa.c)

add_executable(rsa rsa.c)
target_link_libraries(rsa libsa libhgpa)

add_executable(sa1 sa1.c)
target_link_libraries(sa1 libsa libhgpa)

add_executable(sa2 sa2.c)
target_link_libraries(sa2 libsa libhgpa)


add_dependencies(check rsa sa1 sa2)

# Add a test for fms
add_test(NAME sa COMMAND utest.sh ${CMAKE_CURRENT_BINARY_DIR}
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})