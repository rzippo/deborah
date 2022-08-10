cmake_minimum_required(VERSION 3.16)
project(deborah)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "-pedantic-errors")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_CXX_FLAGS "/permissive- /EHsc")
    set(MSVC_COMPILE_FLAGS "/MP")
endif()

add_library(deborah-lib
        src/lib/deborah-globals.h
        src/lib/curve.cpp
        src/lib/curve.h
        src/lib/flow.cpp
        src/lib/flow.h
        src/lib/parametrizedleakybucket.cpp
        src/lib/parametrizedleakybucket.h
        src/lib/linearsegment.cpp
        src/lib/linearsegment.h
        src/lib/nestingtree.cpp
        src/lib/nestingtree.h
        src/lib/node.cpp
        src/lib/node.h
        src/lib/polynomial.cpp
        src/lib/polynomial.h
        src/lib/parametrizedpseudoaffine.cpp
        src/lib/parametrizedpseudoaffine.h
        src/lib/rng.cpp
        src/lib/rng.h
        src/lib/tandem.cpp
        src/lib/tandem.h
        src/lib/timing.h
        src/lib/tnode.cpp
        src/lib/tnode.h
        src/lib/simplex/simplex.cpp
        src/lib/simplex/simplex.h
        src/lib/simplex/constraint.cpp
        src/lib/simplex/constraint.h
        src/lib/simplex/objective.cpp
        src/lib/simplex/objective.h
        src/lib/simplex/solution.cpp
        src/lib/simplex/solution.h
        src/lib/pseudoaffine.cpp
        src/lib/pseudoaffine.h
        src/lib/leakybucket.cpp
        src/lib/leakybucket.h)

add_executable(deborah
        src/deborah-main.cpp
        src/deborah-core.cpp
        src/deborah-core.h
)
target_link_libraries(deborah deborah-lib)

add_executable(deborah-tests
        src/unit-tests/main.cpp
        src/unit-tests/tandem.cpp
        src/unit-tests/simplex.cpp
)
target_link_libraries(deborah-tests deborah-lib)
