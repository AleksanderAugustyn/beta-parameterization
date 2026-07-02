cmake_minimum_required(VERSION 3.20)

include(FetchContent)

project(beta-parameterization
        VERSION 2.0.0
        LANGUAGES Fortran CXX
        DESCRIPTION "Beta parameterization for nuclear physics applications"
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Fetch GCC compiler options package ---
FetchContent_Declare(
        GCC-Compiler-Options
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/GCC-Compiler-Options.git
        GIT_TAG 1.4.1
)
FetchContent_MakeAvailable(GCC-Compiler-Options)

# --- Fetch foundational Fortran modules ---
FetchContent_Declare(
        Fortran-Foundations
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/Fortran-Foundations.git
        GIT_TAG 2.2.2
)
FetchContent_MakeAvailable(Fortran-Foundations)

# --- Preprocessor & build type ---
set(CMAKE_Fortran_PREPROCESS ON)
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "RelWithDebInfo" "Release")
endif ()

# --- Compiler flags via GCC-Compiler-Options ---
include(GCCCompilerOptions/FortranCompilerOptions)
create_fortran_library_interface(TARGET beta_parameterization_flags)

# --- Source files ---
set(BETA_PARAM_SOURCES
        src/beta_parameterization_workers_mod.f08
        src/beta_parameterization_mod.f08
        src/beta_parameterization_c_api_mod.f08
)

# --- Static library target ---
add_library(beta_parameterization STATIC ${BETA_PARAM_SOURCES})
target_link_libraries(beta_parameterization
        PUBLIC fortran_foundations
        PRIVATE beta_parameterization_flags
)
target_compile_definitions(beta_parameterization PRIVATE
        $<$<CONFIG:Debug>:DEBUG>
        $<$<CONFIG:RelWithDebInfo>:RELWITHDEBINFO>
        $<$<CONFIG:Release>:RELEASE>
)
set_target_properties(beta_parameterization PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)
target_include_directories(beta_parameterization
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
add_library(BetaParameterization::beta_parameterization ALIAS beta_parameterization)

# --- Shared library target (for FFI consumers: Python ctypes, Rust, etc.) ---
add_library(beta_parameterization_shared SHARED ${BETA_PARAM_SOURCES})
target_link_libraries(beta_parameterization_shared
        PUBLIC fortran_foundations
        PRIVATE beta_parameterization_flags
)
target_compile_definitions(beta_parameterization_shared PRIVATE
        $<$<CONFIG:Debug>:DEBUG>
        $<$<CONFIG:RelWithDebInfo>:RELWITHDEBINFO>
        $<$<CONFIG:Release>:RELEASE>
)
set_target_properties(beta_parameterization_shared PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files_shared
        OUTPUT_NAME beta_parameterization
)
target_include_directories(beta_parameterization_shared
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files_shared>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
add_library(BetaParameterization::beta_parameterization_shared
        ALIAS beta_parameterization_shared)

# --- C++ header-only wrapper (consumes the shared library) ---
add_library(beta_parameterization_cxx INTERFACE)
target_include_directories(beta_parameterization_cxx INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(beta_parameterization_cxx INTERFACE
        beta_parameterization_shared)
target_compile_features(beta_parameterization_cxx INTERFACE cxx_std_20)
add_library(BetaParameterization::beta_parameterization_cxx
        ALIAS beta_parameterization_cxx)

# --- Tests ---
enable_testing()

add_library(beta_param_test_utils STATIC tests/test_utils_mod.f08)
target_link_libraries(beta_param_test_utils
        PUBLIC fortran_foundations
        PRIVATE beta_parameterization_flags)
set_target_properties(beta_param_test_utils PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files_tests)
target_include_directories(beta_param_test_utils PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files_tests>)

set(BETA_PARAM_TEST_SUITES core com error breakdown golden)
foreach (suite IN LISTS BETA_PARAM_TEST_SUITES)
    add_executable(beta_param_${suite}_test tests/beta_param_${suite}_test.f08)
    target_link_libraries(beta_param_${suite}_test PRIVATE
            beta_parameterization beta_param_test_utils beta_parameterization_flags)
    set_target_properties(beta_param_${suite}_test PROPERTIES
            Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files_tests/${suite})
    add_test(NAME ${suite} COMMAND beta_param_${suite}_test)
endforeach ()

# Golden literals carry 17 significant digits (exact double round-trip);
# gfortran's conversion-extra warning flags the 17th digit as non-significant.
# Source-file property so it lands after the flags-interface options (last wins).
set_source_files_properties(tests/beta_param_golden_test.f08 PROPERTIES
        COMPILE_OPTIONS "-Wno-conversion-extra;-Wno-error=conversion-extra")

add_executable(golden_capture tests/golden_capture.f08)
target_link_libraries(golden_capture PRIVATE
        beta_parameterization beta_param_test_utils beta_parameterization_flags)
set_target_properties(golden_capture PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files_tests/capture)

add_executable(c_api_smoke_test tests/c_api_smoke_test.cpp)
target_link_libraries(c_api_smoke_test PRIVATE beta_parameterization_cxx)
add_test(NAME c_api_smoke COMMAND c_api_smoke_test)

find_package(Threads REQUIRED)
add_executable(thread_stress_test tests/thread_stress_test.cpp)
target_link_libraries(thread_stress_test PRIVATE beta_parameterization_cxx Threads::Threads)
add_test(NAME thread_stress COMMAND thread_stress_test)