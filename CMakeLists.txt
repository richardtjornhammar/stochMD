cmake_minimum_required(VERSION 2.8.9)
project(My_MD_Project)

include_directories(src)
file(GLOB SOURCES "src/*.c")

#SET(GCC_COVERAGE_LINK_FLAGS "-lm")
#SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )

add_executable(stochMD ${SOURCES})
target_link_libraries(stochMD -lm)
