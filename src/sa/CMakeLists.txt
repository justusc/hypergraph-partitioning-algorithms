add_library(ad_lib_sa OBJECT ad_lib_sa.c ad_lib_sa.h)

add_executable(ad_rsa ad_rsa.c $<TARGET_OBJECTS:ad_lib_sa>)
target_link_libraries(ad_rsa ad)

add_executable(ad_sa1 ad_sa1.c $<TARGET_OBJECTS:ad_lib_sa>)
target_link_libraries(ad_sa1 ad)

add_executable(ad_sa2 ad_sa2.c $<TARGET_OBJECTS:ad_lib_sa>)
target_link_libraries(ad_sa2 ad)
