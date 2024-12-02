cmake_minimum_required(VERSION 3.22)
project(complex_predict)

# 启动优化编译
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_STANDARD 17)

# set(CMAKE_C_COMPILER "C:\\APP\\mingw64\\bin\\gcc")
# SET(CMAKE_CXX_COMPILER "C:\\APP\\mingw64\\g++")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# 设置动态库路径
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
# 设置头文件路径
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# 源文件
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
# 库 源文件
file(GLOB LIB_SOURCE ${SOURCE_DIR}/lib/*.cpp)


# 编译动态库
add_library(cpdp SHARED ${LIB_SOURCE})

# # cns
add_executable(pcegs ${SOURCE_DIR}/pcegs.cpp)
target_link_libraries(pcegs cpdp)

## nodynamic
add_executable(nodynamic ${SOURCE_DIR}/no_dynamic.cpp)
target_link_libraries(nodynamic cpdp)

## nodynamic
add_executable(nodag ${SOURCE_DIR}/no_dag.cpp)
target_link_libraries(nodag cpdp)

## nodynamic
add_executable(no_essential ${SOURCE_DIR}/no_essential.cpp)
target_link_libraries(no_essential cpdp)

