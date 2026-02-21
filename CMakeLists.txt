
cmake_minimum_required (VERSION 2.8.11)
project (stateEstimator)

include_directories(inc)

file(GLOB SOURCES "src/*.c")

option(ENABLE_COVERAGE "Enable code coverage instrumentation (GCC/Clang)" OFF)
if(ENABLE_COVERAGE)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        message(STATUS "Coverage enabled (GCC/Clang)")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g --coverage")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    else()
        message(WARNING "ENABLE_COVERAGE is ON, but the compiler does not support GCC/Clang coverage flags.")
    endif()
endif()

option(ENABLE_STRICT_WARNINGS "Enable -Wall/-Wextra and warnings-as-errors" ON)
if(ENABLE_STRICT_WARNINGS)
    if(MSVC)
        add_compile_options(/W4 /WX)
    else()
        add_compile_options(-Wall -Wextra -Werror)
    endif()
endif()

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
add_executable(ekfJacobianUpdateTest test/ekfJacobianUpdateTest.c)
add_executable(ekfAdvancedUpdateTest test/ekfAdvancedUpdateTest.c)
add_executable(ekfFiniteDiffJacobianTest test/ekfFiniteDiffJacobianTest.c)
add_executable(ekfFilterpyExportTest test/ekfFilterpyExportTest.c)
add_executable(ekfFilterpyExportNonlinearTest test/ekfFilterpyExportNonlinearTest.c)
add_executable(ekfMeasurementDimensionTest test/ekfMeasurementDimensionTest.c)
add_executable(ekfPredictUpdateAMatrixNullTest test/ekfPredictUpdateAMatrixNullTest.c)
add_executable(ekfInitAliasAndUpdateUserDataTest test/ekfInitAliasAndUpdateUserDataTest.c)
add_executable(ekfBoatModelTest test/ekfBoatModelTest.c)
add_executable(ekfBoatExportTest test/ekfBoatExportTest.c)
add_executable(kalmanFilterTest test/kalmanFilterTest.c)
add_executable(kalmanFilterExportTest test/kalmanFilterExportTest.c)

if(NOT MSVC)
    set(MATH_LIB m)
endif()

enable_testing()

add_test(NAME ekfMallocTest COMMAND $<TARGET_FILE:ekfMallocTest>)
add_test(NAME ekfStaticInitTest COMMAND $<TARGET_FILE:ekfStaticInitTest>)
add_test(NAME ekfPersonWalkTest COMMAND $<TARGET_FILE:ekfPersonWalkTest>)
add_test(NAME ekfMatrixMathTest COMMAND $<TARGET_FILE:ekfMatrixMathTest>)
add_test(NAME ekfJacobianUpdateTest COMMAND $<TARGET_FILE:ekfJacobianUpdateTest>)
add_test(NAME ekfAdvancedUpdateTest COMMAND $<TARGET_FILE:ekfAdvancedUpdateTest>)
add_test(NAME ekfFiniteDiffJacobianTest COMMAND $<TARGET_FILE:ekfFiniteDiffJacobianTest>)
add_test(NAME ekfFilterpyExportTest COMMAND $<TARGET_FILE:ekfFilterpyExportTest>)
add_test(NAME ekfFilterpyExportNonlinearTest COMMAND $<TARGET_FILE:ekfFilterpyExportNonlinearTest>)
add_test(NAME ekfMeasurementDimensionTest COMMAND $<TARGET_FILE:ekfMeasurementDimensionTest>)
add_test(NAME ekfPredictUpdateAMatrixNullTest COMMAND $<TARGET_FILE:ekfPredictUpdateAMatrixNullTest>)
add_test(NAME ekfInitAliasAndUpdateUserDataTest COMMAND $<TARGET_FILE:ekfInitAliasAndUpdateUserDataTest>)
add_test(NAME ekfBoatModelTest COMMAND $<TARGET_FILE:ekfBoatModelTest>)
add_test(NAME ekfBoatExportTest COMMAND $<TARGET_FILE:ekfBoatExportTest> ${CMAKE_BINARY_DIR}/ekf_boat_compare.csv)
add_test(NAME kalmanFilterTest COMMAND $<TARGET_FILE:kalmanFilterTest>)
add_test(NAME kalmanFilterExportTest COMMAND $<TARGET_FILE:kalmanFilterExportTest> ${CMAKE_BINARY_DIR}/kalman_filterpy_compare.csv)
set_tests_properties(ekfBoatModelTest PROPERTIES WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
set_tests_properties(ekfBoatExportTest PROPERTIES WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
set_tests_properties(kalmanFilterExportTest PROPERTIES WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME ekfMatrixMath_copy COMMAND $<TARGET_FILE:ekfMatrixMathTest> copy)
add_test(NAME ekfMatrixMath_inverse2x2 COMMAND $<TARGET_FILE:ekfMatrixMathTest> inverse2x2)
add_test(NAME ekfMatrixMath_inverse3x3 COMMAND $<TARGET_FILE:ekfMatrixMathTest> inverse3x3)
add_test(NAME ekfMatrixMath_inverse4x4 COMMAND $<TARGET_FILE:ekfMatrixMathTest> inverse4x4)
add_test(NAME ekfMatrixMath_inverse5x5 COMMAND $<TARGET_FILE:ekfMatrixMathTest> inverse5x5)
add_test(NAME ekfMatrixMath_inverse6x6 COMMAND $<TARGET_FILE:ekfMatrixMathTest> inverse6x6)
add_test(NAME ekfMatrixMath_inverse_gauss_jordan COMMAND $<TARGET_FILE:ekfMatrixMathTest> inverse_gauss_jordan)
add_test(NAME ekfMatrixMath_inverse_lu COMMAND $<TARGET_FILE:ekfMatrixMathTest> inverse_lu)
add_test(NAME ekfMatrixMath_inverse_cholesky COMMAND $<TARGET_FILE:ekfMatrixMathTest> inverse_cholesky)
add_test(NAME ekfMatrixMath_inverse_auto_fallback COMMAND $<TARGET_FILE:ekfMatrixMathTest> inverse_auto_fallback)
add_test(NAME ekfMatrixMath_add COMMAND $<TARGET_FILE:ekfMatrixMathTest> add)
add_test(NAME ekfMatrixMath_sub COMMAND $<TARGET_FILE:ekfMatrixMathTest> sub)
add_test(NAME ekfMatrixMath_mult COMMAND $<TARGET_FILE:ekfMatrixMathTest> mult)
add_test(NAME ekfMatrixMath_transpose COMMAND $<TARGET_FILE:ekfMatrixMathTest> transpose)
add_test(NAME ekfMatrixMath_scale COMMAND $<TARGET_FILE:ekfMatrixMathTest> scale)
add_test(NAME ekfMatrixMath_identity_minus_a COMMAND $<TARGET_FILE:ekfMatrixMathTest> identity_minus_a)
add_test(NAME ekfMatrixMath_dimension_mismatch COMMAND $<TARGET_FILE:ekfMatrixMathTest> dimension_mismatch)
add_test(NAME ekfMatrixMath_nan_inf_check COMMAND $<TARGET_FILE:ekfMatrixMathTest> nan_inf_check)
add_test(NAME ekfMatrixMath_null_checks COMMAND $<TARGET_FILE:ekfMatrixMathTest> null_checks)
add_test(NAME ekfMatrixMath_non_init_checks COMMAND $<TARGET_FILE:ekfMatrixMathTest> non_init_checks)
add_test(NAME ekfMatrixMath_alloc_modes COMMAND $<TARGET_FILE:ekfMatrixMathTest> alloc_modes)
add_test(NAME ekfMatrixMath_dim_check_add_macro COMMAND $<TARGET_FILE:ekfMatrixMathTest> dim_check_add_macro)
add_test(NAME ekfMatrixMath_dim_check_mult_macro COMMAND $<TARGET_FILE:ekfMatrixMathTest> dim_check_mult_macro)
add_test(NAME ekfMatrixMath_dim_check_scale_macro COMMAND $<TARGET_FILE:ekfMatrixMathTest> dim_check_scale_macro)
add_test(NAME ekfMatrixMath_static_add COMMAND $<TARGET_FILE:ekfMatrixMathTest> static_add)
add_test(NAME ekfMatrixMath_static_mult COMMAND $<TARGET_FILE:ekfMatrixMathTest> static_mult)
add_test(NAME ekfMatrixMath_static_transpose COMMAND $<TARGET_FILE:ekfMatrixMathTest> static_transpose)
add_test(NAME ekfMatrixMath_static_scale COMMAND $<TARGET_FILE:ekfMatrixMathTest> static_scale)
add_test(NAME ekfMatrixMath_static_identity_minus_a COMMAND $<TARGET_FILE:ekfMatrixMathTest> static_identity_minus_a)
add_test(NAME ekfMatrixMath_compare_mismatch COMMAND $<TARGET_FILE:ekfMatrixMathTest> compare_mismatch)
add_test(NAME ekfMatrixMath_copy_2d_array COMMAND $<TARGET_FILE:ekfMatrixMathTest> copy_2d_array)
add_test(NAME ekfMatrixMath_mult_zeroing COMMAND $<TARGET_FILE:ekfMatrixMathTest> mult_zeroing)
add_test(NAME ekfMatrixMath_inverse_singular COMMAND $<TARGET_FILE:ekfMatrixMathTest> inverse_singular)
add_test(NAME ekfMatrixMath_inverse_jitter COMMAND $<TARGET_FILE:ekfMatrixMathTest> inverse_jitter)
add_test(NAME ekfMatrixMath_nan_propagation COMMAND $<TARGET_FILE:ekfMatrixMathTest> nan_propagation)
add_test(NAME ekfMatrixMath_identity_properties COMMAND $<TARGET_FILE:ekfMatrixMathTest> identity_properties)
add_test(NAME ekfMatrixMath_set_identity COMMAND $<TARGET_FILE:ekfMatrixMathTest> set_identity)
add_test(NAME ekfMatrixMath_transpose_twice COMMAND $<TARGET_FILE:ekfMatrixMathTest> transpose_twice)
add_test(NAME ekfMatrixMath_scale_identity_zero COMMAND $<TARGET_FILE:ekfMatrixMathTest> scale_identity_zero)
add_test(NAME ekfMatrixMath_inverse_sanity COMMAND $<TARGET_FILE:ekfMatrixMathTest> inverse_sanity)
add_test(NAME ekfMatrixMath_inverse_no_mutate COMMAND $<TARGET_FILE:ekfMatrixMathTest> inverse_no_mutate)
add_test(NAME ekfMatrixMath_jagged_static_parity COMMAND $<TARGET_FILE:ekfMatrixMathTest> jagged_static_parity)
add_test(NAME ekfMatrixMath_nonsquare_mult COMMAND $<TARGET_FILE:ekfMatrixMathTest> nonsquare_mult)
add_test(NAME ekfMatrixMath_nan_identity_minus_a COMMAND $<TARGET_FILE:ekfMatrixMathTest> nan_identity_minus_a)
add_test(NAME ekfMatrixMath_aliasing_guards COMMAND $<TARGET_FILE:ekfMatrixMathTest> aliasing_guards)

target_link_libraries(ekfMallocTest stateEstimatorLib ${MATH_LIB})
target_link_libraries(ekfStaticInitTest stateEstimatorLib ${MATH_LIB})
target_link_libraries(ekfPersonWalkTest stateEstimatorLib ${MATH_LIB})
target_link_libraries(ekfMatrixMathTest stateEstimatorLib ${MATH_LIB})
target_link_libraries(ekfJacobianUpdateTest stateEstimatorLib ${MATH_LIB})
target_link_libraries(ekfAdvancedUpdateTest stateEstimatorLib ${MATH_LIB})
target_link_libraries(ekfFiniteDiffJacobianTest stateEstimatorLib ${MATH_LIB})
target_link_libraries(ekfFilterpyExportTest stateEstimatorLib ${MATH_LIB})
target_link_libraries(ekfFilterpyExportNonlinearTest stateEstimatorLib ${MATH_LIB})
target_link_libraries(ekfMeasurementDimensionTest stateEstimatorLib ${MATH_LIB})
target_link_libraries(ekfPredictUpdateAMatrixNullTest stateEstimatorLib ${MATH_LIB})
target_link_libraries(ekfInitAliasAndUpdateUserDataTest stateEstimatorLib ${MATH_LIB})
target_link_libraries(ekfBoatModelTest stateEstimatorLib ${MATH_LIB})
target_link_libraries(ekfBoatExportTest stateEstimatorLib ${MATH_LIB})
target_link_libraries(kalmanFilterTest stateEstimatorLib ${MATH_LIB})
target_link_libraries(kalmanFilterExportTest stateEstimatorLib ${MATH_LIB})

set(TEST_TARGETS
    ekfMallocTest
    ekfStaticInitTest
    ekfPersonWalkTest
    ekfMatrixMathTest
    ekfJacobianUpdateTest
    ekfAdvancedUpdateTest
    ekfFiniteDiffJacobianTest
    ekfFilterpyExportTest
    ekfFilterpyExportNonlinearTest
    ekfMeasurementDimensionTest
    ekfPredictUpdateAMatrixNullTest
    ekfInitAliasAndUpdateUserDataTest
    ekfBoatModelTest
    ekfBoatExportTest
    kalmanFilterTest
    kalmanFilterExportTest
)

option(FORCE_MATRIX_MATH_DEBUG_CHECKS_FOR_TESTS "Force matrix debug checks when running tests (even in Release)" ON)
if(FORCE_MATRIX_MATH_DEBUG_CHECKS_FOR_TESTS)
    # Matrix checks are compiled into the library, so enable them there.
    target_compile_definitions(stateEstimatorLib PRIVATE MATRIX_MATH_DEBUG_CHECKS=1)
    # Also enable for test translation units.
    foreach(tgt IN LISTS TEST_TARGETS)
        target_compile_definitions(${tgt} PRIVATE MATRIX_MATH_DEBUG_CHECKS=1)
    endforeach()
endif()

install(TARGETS ekfMallocTest ekfStaticInitTest ekfPersonWalkTest ekfMatrixMathTest ekfJacobianUpdateTest ekfAdvancedUpdateTest ekfFiniteDiffJacobianTest ekfFilterpyExportTest ekfFilterpyExportNonlinearTest ekfMeasurementDimensionTest ekfPredictUpdateAMatrixNullTest ekfInitAliasAndUpdateUserDataTest ekfBoatModelTest ekfBoatExportTest kalmanFilterTest kalmanFilterExportTest
        DESTINATION bin)

if(ENABLE_COVERAGE AND CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    # Restrict coverage to library sources/headers (exclude tests).
    set(COVERAGE_FILTER_ARGS
        --filter "${CMAKE_SOURCE_DIR}/src/"
        --filter "${CMAKE_SOURCE_DIR}/inc/"
        --exclude "${CMAKE_SOURCE_DIR}/test/"
    )
    add_custom_target(coverage
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
        COMMAND gcovr -r ${CMAKE_SOURCE_DIR} ${COVERAGE_FILTER_ARGS} --html --html-details -o coverage.html
        COMMAND gcovr -r ${CMAKE_SOURCE_DIR} ${COVERAGE_FILTER_ARGS} --xml-pretty -o coverage.xml
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running tests and generating coverage report (gcovr)"
    )
    # Ensure coverage target always builds all test executables before ctest.
    add_dependencies(coverage ${TEST_TARGETS})
endif()
