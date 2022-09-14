cmake_minimum_required(VERSION 3.10)

project(hello)

set(CMAKE_CXX_STANDARD 17)

aux_source_directory(. SRC)

add_executable(main ${SRC})
