#include "camera.h"

#define My_PI (3.141592653589793238462643383279502884)

void camera_t_lockCamera(camera_t* camera)
{
	float camPitch = camera->camPitch;
	if(camPitch>90)
		camPitch=90;
	if(camPitch<-90)
		camPitch=-90;
    camera->camPitch =camPitch;

    float camYaw = camera->camYaw;
	if(camYaw<0.0)
		camYaw+=360.0;
	if(camYaw>360.0)
		camYaw-=360;
    camera->camYaw =camYaw;
}

void camera_t_moveCamera(camera_t* cam, float dir)
{
	float rad=(cam->camYaw+dir)* My_PI /180.0;
	cam->loc.coord[0]-=sin(rad)*cam->movevel;
	cam->loc.coord[2]-=cos(rad)*cam->movevel;
}

void camera_t_moveCameraUp(camera_t* cam, float dir)
{
	float rad=(cam->camPitch+dir)*M_PI/180.0;
	cam->loc.coord[1]+=sin(rad)*cam->movevel;
}

camera_t* camera_t_construct(point_t l)
{
	camera_t* cam = (camera_t*)calloc(1, sizeof(camera_t));
	cam->loc = l;
	cam->camPitch = 0;
	cam->camYaw = 0;
	cam->movevel = 0.2;
	cam->mousevel = 0.2;
	cam->mi = 0;
	cam->ismoved = 0;
	return cam;
}

void camera_t_destruct(camera_t* cam){
    cam->loc = ZERO();
	cam->camPitch = 0;
	cam->camYaw = 0;
	cam->movevel = 0;
	cam->mousevel = 0;
	cam->mi = 0;
	cam->ismoved = 0;
	free(cam);
}

void camera_t_Control(camera_t* cam)
{
	if(cam->mi)
	{
		int MidX=320;
		int MidY=240;
		SDL_ShowCursor(SDL_DISABLE);
		int tmpx,tmpy;
		SDL_GetMouseState(&tmpx,&tmpy);
		cam->camYaw+=cam->mousevel*(MidX-tmpx);
		cam->camPitch+=cam->mousevel*(MidY-tmpy);
		camera_t_lockCamera(cam);
		SDL_WarpMouse(MidX,MidY);
		Uint8* state=SDL_GetKeyState(NULL);
		cam->ismoved=0;
		if(state[SDLK_w])
		{
			cam->ismoved=1;
			if(cam->camPitch!=90 && cam->camPitch!=-90)
				camera_t_moveCamera(cam, 0.0);
			camera_t_moveCameraUp(cam, 0.0);
		}else if(state[SDLK_s])
		{
			cam->ismoved=1;
			if(cam->camPitch!=90 && cam->camPitch!=-90)
				camera_t_moveCamera(cam, 180.0);
			camera_t_moveCameraUp(cam, 180.0);
		}
		if(state[SDLK_a])
		{
			cam->ismoved=1;
			camera_t_moveCamera(cam, 90.0);
		}
		else if(state[SDLK_d])
		{
			cam->ismoved=1;
			camera_t_moveCamera(cam, 270);
		}
	}
	glRotatef(-cam->camPitch,1.0,0.0,0.0);
	glRotatef(-cam->camYaw,0.0,1.0,0.0);
}

void camera_t_UpdateCamera(camera_t* cam)
{
	glTranslatef(-cam->loc.coord[0],-cam->loc.coord[1],-cam->loc.coord[2]);
}

//change the spherical coordinate system to cartesian
point_t camera_t_getVector(camera_t* cam)
{
	//Google->spherical to cartesian
 	return (point_t_get_point(-cos(cam->camPitch*M_PI/180.0)*sin(cam->camYaw*M_PI/180.0),sin(cam->camPitch*M_PI/180.0),-cos(cam->camPitch*M_PI/180.0)*cos(cam->camYaw*M_PI/180.0)));
}

void camera_t_mouseIn(camera_t* cam, int b)
{
	cam->mi=b;
}

int camera_t_isMouseIn(camera_t* cam)
{
	return cam->mi;
}
