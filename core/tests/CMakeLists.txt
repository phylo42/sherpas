cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# Tests
include(CTest)
include(Catch)
enable_testing()

add_executable(tests_n EXCLUDE_FROM_ALL "")

target_sources(tests_n
        PRIVATE
            src/main.cpp
        )

find_package(Boost REQUIRED COMPONENTS filesystem)

target_link_libraries(tests_n
        PRIVATE
            rappas::core_n
            Catch2::Catch2
            Boost::filesystem
        )

catch_discover_tests(tests_n)