cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(rappas_utils LANGUAGES CXX VERSION 0.1.0)

# Create a library target and set properties
add_library(utils STATIC "")
add_library(rappas::utils ALIAS utils)

target_sources(utils
        PRIVATE
        src/io/file_io.cpp
        src/io/fasta.cpp
        src/meta.cpp
        )

find_package(Boost REQUIRED COMPONENTS iostreams)

target_link_libraries(utils
        PUBLIC
            Boost::iostreams
        )

target_include_directories(utils
        PUBLIC
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# Turn on the warnings and treat them as errors
target_compile_options(utils
        PRIVATE
            -Wall -Wextra -Wpedantic -Werror
        )

# explicitly require c++17 standart
# TODO: check this out if it is possible to relax this to PUBLIC cxx_std_14 and PRIVATE cxx_std_17
target_compile_features(utils
        PUBLIC
            cxx_std_17
        )

######################################################################################################
# Library installation
include(GNUInstallDirs)

# Exporting the core target
# 'make install' to the correct locations (provided by GNUInstallDirs)
# TODO: check if it works with 'make install'
install(TARGETS utils
        EXPORT utilsTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

export(TARGETS utils
       FILE ${CMAKE_BINARY_DIR}/rappasConfig.cmake)

export(PACKAGE rappas)