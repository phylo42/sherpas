cmake_minimum_required(VERSION 3.10)

add_executable(
	rappas-recombn
        src/main.cpp
        src/arcs.h
        src/arcs.cpp
        src/heap.h
        src/heap.cpp
        src/query.h
        src/query.cpp
        src/write.cpp
        src/write.h
)

find_package(Boost REQUIRED COMPONENTS program_options filesystem iostreams)
find_package(OpenMP REQUIRED)

target_link_libraries(rappas-recombn
        PRIVATE
            rappas::core_n
            Boost::iostreams
            OpenMP::OpenMP_CXX
        )
