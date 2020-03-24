cmake_minimum_required(VERSION 3.5)

project(mvp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# change accordingly
include_directories("/usr/local/include")

# nlopt has to be installed, change accordingly
find_library(LNLOPT libnlopt.so "/usr/local/lib" )

add_executable(covid covid.cpp)

target_link_libraries("covid" ${LNLOPT} )


