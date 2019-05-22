cmake_minimum_required(VERSION 3.10)
project(rappas-recomb)

# default build type Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# define a hash map
if(NOT HASH_MAP)
    set(HASH_MAP "USE_TSL_HOPSCOTCH_MAP")
endif()

add_subdirectory(core)
add_subdirectory(recomb)
