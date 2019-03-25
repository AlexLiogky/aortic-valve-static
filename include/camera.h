#ifndef CAMERA_H_INCLUDED
#define CAMERA_H_INCLUDED

#include <math.h>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include "point.h"


typedef struct camera_t{
    point_t loc;
	float camPitch,camYaw;
	float movevel;
	float mousevel;
	int mi,ismoved;
} camera_t;

void camera_t_lockCamera(camera_t* camera);

void camera_t_moveCamera(camera_t* cam, float dir);

void camera_t_moveCameraUp(camera_t* cam, float dir);

camera_t* camera_t_construct(point_t l);

void camera_t_destruct(camera_t* cam);

void camera_t_Control(camera_t* cam);

void camera_t_UpdateCamera(camera_t* cam);

//change the spherical coordinate system to cartesian
point_t camera_t_getVector(camera_t* cam);

void camera_t_mouseIn(camera_t* cam, int b);

int camera_t_isMouseIn(camera_t* cam);

#endif // CAMERA_H_INCLUDED
