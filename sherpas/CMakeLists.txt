cmake_minimum_required(VERSION 3.10)

add_executable(
	SHERPAS
        src/main.cpp
        src/arcs.h
        src/arcs.cpp
        src/heap.h
        src/heap.cpp
        src/query.h
        src/query.cpp
        src/groups.cpp
        src/groups.h
        src/output.cpp
        src/output.h
)

target_link_libraries(SHERPAS
        PRIVATE
            xpas::dna
        )
