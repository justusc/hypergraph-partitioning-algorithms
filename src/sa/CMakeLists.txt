add_library(libsa ad_lib_sa.c)

add_executable(ad_rsa ad_rsa.c)
target_link_libraries(ad_rsa libsa libhgpa)

add_executable(ad_sa1 ad_sa1.c)
target_link_libraries(ad_sa1 libsa libhgpa)

add_executable(ad_sa2 ad_sa2.c)
target_link_libraries(ad_sa2 libsa libhgpa)


add_dependencies(check ad_rsa ad_sa1 ad_sa2)

# Add a test for fms
add_test(NAME sa COMMAND utest.sh ${CMAKE_CURRENT_BINARY_DIR}
         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})