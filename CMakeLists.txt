cmake_minimum_required(VERSION 3.14)

project(DPG)

set(CMAKE_CXX_STANDARD 17)

# Finding the pz package
find_package(NeoPZ REQUIRED HINTS ../NeoPZ_install ../neopz_install)

# Mac OS X specific code
if (APPLE)
  add_definitions(-DMACOSX)
endif()
# ----------------

add_library(DPG_lib SHARED "")
include (cmake/EnableEigen.cmake)
enable_eigen(DPG_lib)
target_link_libraries(DPG_lib PUBLIC NeoPZ::pz)
target_include_directories(DPG_lib PUBLIC ${PZ_INCLUDE_DIRS})
set_target_properties(DPG_lib PROPERTIES LINKER_LANGUAGE CXX)

add_subdirectory(src)
# get_target_property(lib_files DPG_lib SOURCES)
# source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${lib_files})

add_subdirectory(targets)