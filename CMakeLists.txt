cmake_minimum_required(VERSION 3.8...3.15)

project (JC_Code LANGUAGES CXX)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# this option is useful for static analyzers
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# Finding the pz package if not in neopz buildtree
if(NOT CMAKE_IS_PZ_BUILDTREE)
    find_package(NeoPZ REQUIRED)
endif()

add_subdirectory(LCCCommon)
add_subdirectory(LCCAdaptive)
add_subdirectory(LCCEstimator)

include_directories(LCCCommon)
include_directories(LCCAdaptive)
include_directories(LCCEstimator)

add_subdirectory(LCC_Projects)

