
cmake_minimum_required (VERSION 2.8.11)
project (stateEstimator)

include_directories(inc)

file(GLOB SOURCES "src/*.c")

add_library(stateEstimatorLib STATIC ${SOURCES})
set_target_properties(stateEstimatorLib PROPERTIES OUTPUT_NAME stateEstimatorLib PUBLIC_HEADER inc/EKF.h)

install(TARGETS stateEstimatorLib 
        DESTINATION ${CMAKE_SOURCE_DIR}/lib)

#testing binaries
add_executable(ekfMallocTest test/ekfMallocTest.c)
add_executable(ekfStaticInitTest test/ekfStaticInitTest.c)
add_executable(ekfPersonWalkTest test/ekfPersonWalkTest.c)

enable_testing()

add_test(NAME ekfMallocTest COMMAND $<TARGET_FILE:ekfMallocTest>)
add_test(NAME ekfStaticInitTest COMMAND $<TARGET_FILE:ekfStaticInitTest>)

target_link_libraries(ekfMallocTest stateEstimatorLib)
target_link_libraries(ekfStaticInitTest stateEstimatorLib)
target_link_libraries(ekfPersonWalkTest stateEstimatorLib)

install(TARGETS ekfMallocTest ekfStaticInitTest ekfPersonWalkTest
        DESTINATION bin)