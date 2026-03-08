cmake_minimum_required(VERSION 3.20)

include(FetchContent)

project(beta_parameterization
        VERSION 1.0.0
        LANGUAGES Fortran
        DESCRIPTION "Beta parameterization for nuclear physics applications"
)

# --- Fetch GCC compiler options package ---
FetchContent_Declare(
        GCC-Compiler-Options
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/GCC-Compiler-Options.git
        GIT_TAG 1.3.0
)
FetchContent_MakeAvailable(GCC-Compiler-Options)

# --- Fetch foundational Fortran modules ---
FetchContent_Declare(
        Fortran-Foundations
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/Fortran-Foundations.git
        GIT_TAG 2.1.0
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

# --- Library target ---
add_library(beta_parameterization STATIC src/beta_parameterization_mod.f08)
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
        PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files>
)
add_library(BetaParameterization::beta_parameterization ALIAS beta_parameterization)