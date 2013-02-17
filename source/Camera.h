// Contributors Justin Hutchison yibbidy@gmail.com

#ifndef Camera_h
#define Camera_h

#include <math.h>
#include <gl/gl.h>
#include "gl/glut.h"
#include "Geometry.h"

class Camera 
{
public:
	double fTheta;
	double fPhi;
	float fRadius; // larger means further away
	bool fRightButtonDown;
	bool fLeftButtonDown;
	int fLastX;
	int fLastY;

	Camera() 
	{
		fTheta = 3.141592654/4.0;
		fPhi = 3.141592654;
		fRadius = 1;
		fRightButtonDown = false;
		fLeftButtonDown = false;
		fLastX = 0;
		fLastY = 0;
	}

	TVector<float, 3> Location() 
	{
		return TVector<float, 3>((float) (-sin(fTheta)*sin(fPhi)), (float) cos(fTheta), (float) (-sin(fTheta)*cos(fPhi)))*fRadius;
	}
	
	TVector<float, 3> Direction() 
	{
		return -TVector<float, 3>((float) (-sin(fTheta)*sin(fPhi)), (float) cos(fTheta), (float) (-sin(fTheta)*cos(fPhi)));
	}
	
	void gluLookAt() 
	{
		float camera[3];
		camera[0] = (float) (-sin(fTheta)*sin(fPhi));
		camera[2] = (float) (-sin(fTheta)*cos(fPhi));
		camera[1] = (float) cos(fTheta);

		float up[3];
		up[0] = (float) (cos(fTheta)*sin(fPhi));
		up[2] = (float) (cos(fTheta)*cos(fPhi));
		up[1] = (float) sin(fTheta);

		::gluLookAt(fRadius*camera[0], fRadius*camera[1], fRadius*camera[2], 0, 0, 0, up[0], up[1], up[2]);
	}


	void Display() 
	{
		GLboolean lighting = glIsEnabled(GL_LIGHTING);
		if( lighting ) glDisable(GL_LIGHTING);

		glBegin(GL_LINES);
		glColor3f(1, 0, 0);
		glVertex3f(0, 0, 0);
		glVertex3f(1, 0, 0);
	
		glColor3f(0, 1, 0);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 1, 0);

		glColor3f(0, 0, 1);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, 1);
		glEnd();

		if( lighting ) glEnable(GL_LIGHTING);
	}


	void MouseFunc(int button, int state, int x, int y) 
	{
		if( GLUT_DOWN == state ) 
		{
			fLastX = x;
			fLastY = y;

			if( button == GLUT_LEFT_BUTTON ) fLeftButtonDown = true;
			else if( button == GLUT_RIGHT_BUTTON ) fRightButtonDown = true;
		} 
		else if( GLUT_UP == state ) 
		{
			if( button == GLUT_LEFT_BUTTON ) fLeftButtonDown = false;
			else if( button == GLUT_RIGHT_BUTTON ) fRightButtonDown = false;
		}


	}

	void MotionFunc(int x, int y)
	{
		if( fLeftButtonDown ) 
		{
			fTheta += .01 * (fLastY-y);
			fPhi += .01 * (fLastX-x);
		}

		if( fRightButtonDown ) 
		{
			fRadius += .01 * (y-fLastY);
		}
		
		fLastX = x;
		fLastY = y;

		glutPostRedisplay();
	}

};


#endif
