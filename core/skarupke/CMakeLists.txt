cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(skarupke VERSION 1.0.0 LANGUAGES CXX)

add_library(flat_hash_map INTERFACE)
add_library(skarupke::flat_hash_map ALIAS flat_hash_map)

target_include_directories(flat_hash_map
        INTERFACE
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        )

###################################################################################################
# Installation instructions

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/skarupke)

install(TARGETS flat_hash_map
        EXPORT skarupke-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Export the target
install(EXPORT skarupke-targets
        FILE skarupkeTargets.cmake
        NAMESPACE skarupke::
        DESTINATION ${INSTALL_CONFIGDIR}
        )

export(EXPORT skarupke-targets
        FILE ${CMAKE_CURRENT_BINARY_DIR}/skarupkeTargets.cmake
        NAMESPACE skarupke::)

export(PACKAGE skarupke)