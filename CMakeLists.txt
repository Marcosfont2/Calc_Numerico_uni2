cmake_minimum_required(VERSION 3.10)
project(TrabalhoCalculoNumerico)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

add_executable(trabalho
    src/main.cpp
    src/gauss.cpp
    src/lu.cpp
    src/jacobi.cpp
    src/operations.cpp
    src/io.cpp
)

add_executable(teste_interpolacao tests/teste_interpolacao.cpp src/gauss.cpp src/operations.cpp)
