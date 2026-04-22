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

# --- C++ header-only wrapper (consumes the shared library) ---
add_library(beta_parameterization_cxx INTERFACE)
target_include_directories(beta_parameterization_cxx INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(beta_parameterization_cxx INTERFACE
        beta_parameterization_shared)
target_compile_features(beta_parameterization_cxx INTERFACE cxx_std_20)
add_library(BetaParameterization::beta_parameterization_cxx
        ALIAS beta_parameterization_cxx)