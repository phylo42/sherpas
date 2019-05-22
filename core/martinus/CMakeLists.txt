cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(martinus VERSION 1.0.0 LANGUAGES CXX)

add_library(robin_hood_map INTERFACE)
add_library(martinus::robin_hood_map ALIAS robin_hood_map)

target_include_directories(robin_hood_map
        INTERFACE
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/robin-hood-hashing/src/include>
        )

###################################################################################################
# Installation instructions

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/martinus)

install(TARGETS robin_hood_map
        EXPORT martinus-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )

install(DIRECTORY robin-hood-hashing/src/include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Export the target
install(EXPORT martinus-targets
        FILE martinusTargets.cmake
        NAMESPACE martinus::
        DESTINATION ${INSTALL_CONFIGDIR}
        )

export(EXPORT martinus-targets
        FILE ${CMAKE_CURRENT_BINARY_DIR}/martinusTargets.cmake
        NAMESPACE martinus::)

export(PACKAGE martinus)