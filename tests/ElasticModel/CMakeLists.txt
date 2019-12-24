cmake_minimum_required(VERSION 2.6)

project(ElastModelTest)
enable_language (CXX C Fortran)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_C_COMPILER gcc)
set(CMAKE_Fortran_COMPILER gfortran)

set (ANI_DIR "/home/alex/Desktop/Work/ani3D-3.1")
set (BULLET_DIR "/home/alex/Desktop/prog/C++-prog/Bullet/BulletTrainig/bullet3-2.88")

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS} -O0 -g")
set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_FORTRAN_FLAGS_DEBUG} -O0 -g")
set(CMAKE_Fortran_FLAGS_RELEASE  "${CMAKE_FORTRAN_FLAGS_RELEASE} -O3")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0  -g")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0  -g")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0  -g -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0  -g -std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -std=c++11")


link_directories(${ANI_DIR}/lib)
link_directories(${BULLET_DIR}/bin)
#include_directories(${ANI_DIR}/include)
include_directories(${ANI_DIR}/src/aniAFT/src/aniFRT/PRM)
include_directories(../../include)
include_directories(${BULLET_DIR}/src)




set(EXEC_NAME model-test)
set(SRC "../../src")
#list(APPEND SOURCE_FILES ${SRC}/ani_exp_stress.c)
#list(APPEND SOURCE_FILES ${SRC}/ani_lin_stress.c)
#list(APPEND SOURCE_FILES ${SRC}/anisotrop-stress.c)
list(APPEND SOURCE_FILES ${SRC}/bound-box.c)
list(APPEND SOURCE_FILES ${SRC}/camera.c)
list(APPEND SOURCE_FILES ${SRC}/computation.c)
list(APPEND SOURCE_FILES ${SRC}/data.c)
list(APPEND SOURCE_FILES ${SRC}/elem.c)
list(APPEND SOURCE_FILES ${SRC}/format_in.c)
list(APPEND SOURCE_FILES ${SRC}/hermit-splines.c)
list(APPEND SOURCE_FILES ${SRC}/intersection.c)
list(APPEND SOURCE_FILES ${SRC}/isotrop-stress.c)
list(APPEND SOURCE_FILES ${SRC}/nets.c)
list(APPEND SOURCE_FILES ${SRC}/node.c)
list(APPEND SOURCE_FILES ${SRC}/point.c)
list(APPEND SOURCE_FILES ${SRC}/postcomputation.c)
list(APPEND SOURCE_FILES ${SRC}/precomputation.c)
list(APPEND SOURCE_FILES ${SRC}/save-data.c)
list(APPEND SOURCE_FILES ${SRC}/separate.c)
list(APPEND SOURCE_FILES ${SRC}/sew_leaf.c)
list(APPEND SOURCE_FILES ${SRC}/spring.c)

#list(APPEND SOURCE_FILES ${SRC}/NetCollisionShape.cpp)
#list(APPEND SOURCE_FILES ${SRC}/NetConcaveCollisionAlgorithm.cpp)
#list(APPEND SOURCE_FILES ${SRC}/NetNetCollisionAlgorithm.cpp)
#list(APPEND SOURCE_FILES ${SRC}/Net_Wraper.cpp)
#list(APPEND SOURCE_FILES ${SRC}/Solver.cpp)
#list(APPEND SOURCE_FILES ${SRC}/StaticCollisionAlgorithm.cpp)
#list(APPEND SOURCE_FILES ${SRC}/World.cpp)

list(APPEND SOURCE_FILES main.cpp)
add_executable(${EXEC_NAME} ${SOURCE_FILES} )

target_link_libraries(${EXEC_NAME} -lfrtprm-3.0 -lmba3D-3.0 -laft3D-3.0 -llapack-3.2 -llapack_ext-3.2 -lblas -lgfortran -lm)
target_link_libraries(${EXEC_NAME} -lSDL -lGL -lGLU -lSDL_image)
target_link_libraries(${EXEC_NAME} -lgsl -lgslcblas)

