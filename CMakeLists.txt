
cmake_minimum_required (VERSION 2.8.11)
project (stateEstimator)

include_directories(inc)

file(GLOB SOURCES "src/*.c")

#set(CMAKE_C_FLAGS "-Wall -Wextra -Werror -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function -Wno-unused-value -Wno-unused-result -Wno-unused-label -Wno-unused-local-typedefs -Wno-unused-macros -Wno-u")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__FILENAME__='\"$(subst
#  ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

add_library(stateEstimatorLib STATIC ${SOURCES})
set_target_properties(stateEstimatorLib PROPERTIES OUTPUT_NAME stateEstimatorLib PUBLIC_HEADER inc/EKF.h)

install(TARGETS stateEstimatorLib 
        DESTINATION ${CMAKE_SOURCE_DIR}/lib)

#testing binaries
add_executable(ekfMallocTest test/ekfMallocTest.c)
add_executable(ekfStaticInitTest test/ekfStaticInitTest.c)
add_executable(ekfPersonWalkTest test/ekfPersonWalkTest.c)
add_executable(ekfMatrixMathTest test/ekfMatrixMathTest.c)

enable_testing()

add_test(NAME ekfMallocTest COMMAND $<TARGET_FILE:ekfMallocTest>)
add_test(NAME ekfStaticInitTest COMMAND $<TARGET_FILE:ekfStaticInitTest>)
add_test(NAME ekfPersonWalkTest COMMAND $<TARGET_FILE:ekfPersonWalkTest>)
add_test(NAME ekfMatrixMathTest COMMAND $<TARGET_FILE:ekfMatrixMathTest>)

target_link_libraries(ekfMallocTest stateEstimatorLib)
target_link_libraries(ekfStaticInitTest stateEstimatorLib)
target_link_libraries(ekfPersonWalkTest stateEstimatorLib)
target_link_libraries(ekfMatrixMathTest stateEstimatorLib)

install(TARGETS ekfMallocTest ekfStaticInitTest ekfPersonWalkTest ekfMatrixMathTest
        DESTINATION bin)