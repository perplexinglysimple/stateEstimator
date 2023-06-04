
cmake_minimum_required (VERSION 2.8.11)
project (stateEstimator)

include_directories(inc)

file(GLOB SOURCES "src/*.c")

file(GLOB TESTS "test/*.c")

add_library(stateEstimatorLib SHARED ${SOURCES})

add_executable(stateEstimatorTests ${TESTS})

install(TARGETS stateEstimatorLib DESTINATION lib)