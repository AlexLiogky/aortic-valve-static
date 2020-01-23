cmake_minimum_required(VERSION 2.6)

project(Aortic-Valve)
enable_language (CXX C Fortran)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_C_COMPILER gcc)
set(CMAKE_Fortran_COMPILER gfortran)

set (ANI_DIR "/home/alex/Desktop/Work/ani3D-3.1")
set (BULLET_DIR "/home/alex/Desktop/prog/C++-prog/Bullet/BulletTrainig/bullet3-2.88")

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS} -O0 -g")
set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_FORTRAN_FLAGS_DEBUG} -O0 -g")
set(CMAKE_Fortran_FLAGS_RELEASE  "${CMAKE_FORTRAN_FLAGS_RELEASE} -O3")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0  -g")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0  -g")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0  -g -std=c++17")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0  -g -std=c++17")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -std=c++17")

#message("build_type = ${CMAKE_BUILD_TYPE}")


link_directories(${ANI_DIR}/lib)
link_directories(${BULLET_DIR}/bin)
#include_directories(${ANI_DIR}/include)
include_directories(${ANI_DIR}/src/aniAFT/src/aniFRT/PRM)
include_directories(include)
include_directories(include/elastic_forces)
#include_directories(libs/LEDA/incl)
#include_directories(include/eigen include/autodiff include/MinEnergyDeformator.h)
include_directories(${BULLET_DIR}/src)




set(EXEC_NAME Aortic-Valve)
set(SRC "src")
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
list(APPEND SOURCE_FILES ${SRC}/elasic_forces/MSM.c)
list(APPEND SOURCE_FILES ${SRC}/elasic_forces/TRQS.c)
list(APPEND SOURCE_FILES ${SRC}/elasic_forces/NEOGOOK.c)
list(APPEND SOURCE_FILES ${SRC}/elasic_forces/REINFORCING.c)
list(APPEND SOURCE_FILES ${SRC}/elasic_forces/elastic_forces.c)

list(APPEND SOURCE_FILES ${SRC}/NetCollisionShape.cpp)
list(APPEND SOURCE_FILES ${SRC}/NetConcaveCollisionAlgorithm.cpp)
list(APPEND SOURCE_FILES ${SRC}/NetNetCollisionAlgorithm.cpp)
list(APPEND SOURCE_FILES ${SRC}/Net_Wraper.cpp)
#list(APPEND SOURCE_FILES ${SRC}/Solver.cpp)
list(APPEND SOURCE_FILES ${SRC}/StaticCollisionAlgorithm.cpp)
list(APPEND SOURCE_FILES ${SRC}/World.cpp)
list(APPEND SOURCE_FILES ${SRC}/MinEnergyDeformator.cpp)
list(APPEND SOURCE_FILES ${SRC}/InputProcessor.cpp)
list(APPEND SOURCE_FILES ${SRC}/Meassure.cpp)

list(APPEND SOURCE_FILES ${SRC}/main.cpp)
add_executable(${EXEC_NAME} ${SOURCE_FILES} )

#set(BUG_NAME buger)
#add_executable(${BUG_NAME} test.cpp)
target_link_libraries(${EXEC_NAME} -lfrtprm-3.0 -lmba3D-3.0 -laft3D-3.0 -llapack-3.2 -llapack_ext-3.2 -lblas -lgfortran -lm -lstdc++fs)
target_link_libraries(${EXEC_NAME} -lSDL -lGL -lGLU -lSDL_image)
target_link_libraries(${EXEC_NAME} -lgsl -lgslcblas)
target_link_libraries(${EXEC_NAME} -lgmpxx -lgmp -lmpfr -lCGAL)
#target_link_libraries(${EXEC_NAME} /home/alex/Desktop/MV/valve model-static/Aortic-Valve/libs/LEDA/libleda.a)
if ( CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "debug")
	target_link_libraries(${EXEC_NAME} -lBulletCollision_codeblocks_x64_debug -lLinearMath_codeblocks_x64_debug -lBullet3Collision_codeblocks_x64_debug)
else()
	target_link_libraries(${EXEC_NAME} -lBulletCollision_codeblocks_x64_release -lLinearMath_codeblocks_x64_release -lBullet3Collision_codeblocks_x64_release)
endif()





