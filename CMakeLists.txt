cmake_minimum_required(VERSION 3.20)

include(FetchContent)

project(beta-parameterization
        VERSION 4.0.0
        LANGUAGES Fortran CXX
        DESCRIPTION "Beta parameterization for nuclear physics applications"
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Fetch GCC compiler options package ---
FetchContent_Declare(
        gcc-compiler-options
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/gcc-compiler-options.git
        GIT_TAG 2.0.0
)
FetchContent_MakeAvailable(gcc-compiler-options)

# --- Fetch foundational Fortran modules ---
FetchContent_Declare(
        fortran-foundations
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/fortran-foundations.git
        GIT_TAG 3.0.0
)
FetchContent_MakeAvailable(fortran-foundations)

# --- Preprocessor & build type ---
set(CMAKE_Fortran_PREPROCESS ON)
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "RelWithDebInfo" "Release")
endif ()

# --- Compiler flags via gcc-compiler-options ---
include(GCCCompilerOptions/FortranCompilerOptions)
create_fortran_library_interface(TARGET beta_parameterization_flags)

# --- Source files ---
set(BETA_PARAM_SOURCES
        src/beta_parameterization_workers_mod.f08
        src/beta_parameterization_mod.f08
        src/beta_parameterization_c_api_mod.f08
)

# --- No LTO for the library's own objects ---
# Release builds use -flto=auto -ffast-math. With LTO the summation kernels are
# inlined into every caller and optimized per call site, so the same call can
# return different last bits from different call sites — which breaks the
# contract's bitwise guarantees (one-shot == cached, statelessness,
# concurrency). Compiled without LTO, each kernel is ONE machine-code body that
# every caller executes. Consumers keep LTO for their own code.
# Source-file property so it lands after the flags-interface options (last wins).
set_source_files_properties(${BETA_PARAM_SOURCES} PROPERTIES COMPILE_OPTIONS "-fno-lto")

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

# --- Install the shared library into the Python package (wheel builds) ---
install(TARGETS beta_parameterization_shared
        LIBRARY DESTINATION beta_parameterization)

# --- Tests (skipped for wheel builds via -DBETA_PARAM_BUILD_TESTS=OFF) ---
option(BETA_PARAM_BUILD_TESTS "Build the Fortran/C++ test suite" ON)
if (BETA_PARAM_BUILD_TESTS)
    enable_testing()

add_library(beta_param_test_utils STATIC tests/test_utils_mod.f08)
target_link_libraries(beta_param_test_utils
        PUBLIC fortran_foundations
        PRIVATE beta_parameterization_flags)
set_target_properties(beta_param_test_utils PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files_tests)
target_include_directories(beta_param_test_utils PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files_tests>)

set(BETA_PARAM_TEST_SUITES status lifecycle resolve outputs standalone property golden
        equivalence statelessness boundary)
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

# Capture tool for the golden baseline: an executable, never a test.
add_executable(golden_capture tests/golden_capture.f08)
target_link_libraries(golden_capture PRIVATE
        beta_parameterization beta_param_test_utils beta_parameterization_flags)
set_target_properties(golden_capture PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files_tests/capture)

# PES-range sweep: 1.55M shapes, far too slow for ctest. An executable driven
# by an explicit custom target, never registered with add_test.
add_executable(beta_pes_sweep tests/beta_pes_sweep.f08)
target_link_libraries(beta_pes_sweep PRIVATE
        beta_parameterization beta_param_test_utils beta_parameterization_flags)
set_target_properties(beta_pes_sweep PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files_tests/sweep)
add_custom_target(test_beta_pes_sweep
        COMMAND beta_pes_sweep
        DEPENDS beta_pes_sweep
        COMMENT "Running the PES-range sweep (1.55M + 172.8k shapes)")

# C API smoke test: links the SHARED library, the binary the bindings load.
add_executable(c_api_smoke_test tests/c_api_smoke_test.cpp)
target_link_libraries(c_api_smoke_test PRIVATE beta_parameterization_cxx)
add_test(NAME c_api_smoke COMMAND c_api_smoke_test)

# Contract family 3 (concurrency): many threads on one shared cache, via the
# C++ wrapper.
find_package(Threads REQUIRED)
add_executable(thread_stress_test tests/thread_stress_test.cpp)
target_link_libraries(thread_stress_test PRIVATE
        beta_parameterization_cxx Threads::Threads)
add_test(NAME thread_stress COMMAND thread_stress_test)

endif ()