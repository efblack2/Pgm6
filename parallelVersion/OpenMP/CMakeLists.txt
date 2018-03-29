cmake_minimum_required(VERSION 3.9)

## set appropriately the following line to create a Double Precision version
set(DOUBLE_PRECISION ON)

project(streamOpenMP)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)


## the following lines mainly contain flags I have not found in CMake
if (CMAKE_C_COMPILER_ID STREQUAL "GNU" )
    #MESSAGE( STATUS "                     yes, the compiler is: "            ${CMAKE_C_COMPILER_ID} )
    set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Wall -Wextra  -ftree-vectorize -fopt-info-vec-optimized -ffast-math -mavx2 ")
    set(CMAKE_EXE_LINKER_FLAGS "-lm" )
elseif ( CMAKE_C_COMPILER_ID STREQUAL "Intel" )
    set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Wall -qopt-report=3 -qopt-report-phase=vec -xCORE-AVX2")
elseif ( CMAKE_C_COMPILER_ID STREQUAL "PGI" )
    set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Minform=inform -Minfo=all -Mvect=simd ")
endif ( )
## end of the following lines mainly contain flags I have not found in CMake

if (DOUBLE_PRECISION)
    set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -DDOUBLE")
endif()


find_package(OpenMP REQUIRED)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")


include_directories("src/include")
add_subdirectory(src)

add_executable(p6 ${p6_Src})
target_link_libraries(p6 ${CMAKE_EXE_LINKER_FLAGS})

############### Debug messages ####################
#MESSAGE( STATUS "CMAKE_BINARY_DIR: "            ${CMAKE_BINARY_DIR} )
#MESSAGE( STATUS "CMAKE_C_COMPILER_ID: "         ${CMAKE_C_COMPILER_ID} )
#MESSAGE( STATUS "CMAKE_C_FLAGS: "               ${CMAKE_C_FLAGS} )
#MESSAGE( STATUS "CMAKE_CURRENT_SOURCE_DIR: "    ${CMAKE_CURRENT_SOURCE_DIR} )
########### end of Debug messages #################

