cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

# Create an executable target and set properties
add_executable(core_example_search EXCLUDE_FROM_ALL "")

target_sources(core_example_search
        PRIVATE
            src/main.cpp
        )

target_link_libraries(core_example_search
        PRIVATE
            rappas::core_n
        )

# Turn on the warnings and treat them as errors
target_compile_options(core_example_search
        PRIVATE
            -Wall -Wextra -Wpedantic -Werror
        )

target_compile_features(core_example_search
        PRIVATE
            cxx_std_17
        )