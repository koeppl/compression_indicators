cmake_minimum_required(VERSION 3.0.2)
project (compression_indicators)


if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

set(CXX_STANDARD gnu++17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -march=native -std=${CXX_STANDARD} ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mtune=native")

include_directories(external/separate_chaining/include)
include_directories(external/bit_span/include)
link_libraries(stdc++fs)

##########
# executables
##########

add_executable (effective_alphabet effective_alphabet.cpp)
add_executable (empirical_entropy empirical_entropy.cpp)
add_executable (h0 h0.cpp)
add_executable (h1 h1.cpp)
add_executable (h2 h2.cpp)
# target_link_libraries(empirical_entropy stdc++fs)
# target_link_libraries(h0 stdc++fs)
# target_link_libraries(h1 stdc++fs)
# target_link_libraries(h2 stdc++fs)


MESSAGE( STATUS "Built Type: " ${CMAKE_BUILD_TYPE} )
MESSAGE( STATUS "CXX Flags: " ${CMAKE_CXX_FLAGS} )
