// Contributors Justin Hutchison yibbidy@gmail.com

#include <stdlib.h>
#include "gl/glut.h"
#include <windows.h>
#include <math.h>
#include "Camera.h"
#include <vector>
#include <strstream>
#include <iomanip>

int gWindowWidth;
int gWindowHeight;

Camera gCamera;

double gOriginalAngle = 45;
double gAngle = 45;

int gBoxCount = 3;

int gAnimationState = 0;
const int gRandomize = 0;
const int gDeRotate = 1;
const int gComputeBounds = 2;
const int gRotateBack = 3;
const int gDone = 4;

int gCurrentBox;
int gBoxEdge; // 0==left 1==bottom 2==right 3==top
				


int round(double d) {
	return (d<.5) ? (int)d : ((int)d) + 1;
}

int RandomInt(int min, int max) {
	return round(rand()/(double)(RAND_MAX)*(max-min)) + min;
}

void InitGraphics(void) {
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void Draw4x4Matrix(float * matrix, float initY) {
	glColor3f(1, 1, 1);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, gWindowWidth, 0, gWindowHeight, -1, 1);
	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	glBegin(GL_QUADS);
	glColor4f(0, 0, 0, .4);
	glVertex2f(0, initY+15);
	glVertex2f(150+70*4, initY+15);
	glVertex2f(150+70*4, initY+15-15*4);
	glVertex2f(0, initY+15-15*4);
	glEnd();

	glDisable(GL_BLEND);


	{
		glColor3f(1, 1, 1);
		glRasterPos2i(0, initY-15);

		std::string text = "ViewMatrix = ";
		for( int c=0; c<(int)text.size(); ++c ) {
			glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, text[c]);
		}
		
	}

	int x = 150;

	for( int i=0; i<4; ++i ) {
		int y = initY;

		if( 0 == i ) glColor3f(1, .2, .2);
		else if( 1 == i ) glColor3f(.2, 1, .2);
		else if( 2 == i ) glColor3f(.2, .2, 1);
		else glColor3f(1, 1, 1);

		for( int j=0; j<4; ++j ) {
			if( 3 == j ) glColor3f(1, 1, 1);

			std::ostrstream str;
			double element = *matrix++;
			if( abs(element) < 1e-2 ) element = 0;
			str << std::setfill('0') << std::setprecision(3) << element << '\0';

			std::string text = str.str();
			
			glRasterPos2i(x, y);
			for( int c=0; c<(int)text.size(); ++c ) {
				glutBitmapCharacter(GLUT_BITMAP_9_BY_15, text[c]);
			}
			glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ' ');
			y -= 15;

		}

		x += 70;
	}
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);

}

typedef TVector<double, 2> Vec;
struct Box {
	Box() {
	}

	Box(Vec bottomLeft, double width, double height) {
		this->bottomLeft = bottomLeft;
		this->width = width;
		this->height = height;
	}

	Vec bottomLeft;
	double width, height;

	float r, g, b;
};

Box gBoxes[13];
Box gBounds;

void InitData() {
	gBoxes[0] = Box(Vec(0, -20), 50, 70);
	gBoxes[1] = Box(Vec(-30, 30), 20, 20);
	gBoxes[2] = Box(Vec(40, 60), 40, 35);

}


void DisplayFunc(void) {
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();

	glRotatef(gAngle, 0, 0, 1);

	for( int i=0; i<gBoxCount; ++i ) {
		Box & box = gBoxes[i];

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		
		glColor3f(box.r*.5, box.g*.5, box.b*.5);
		glBegin(GL_QUADS);
		glVertex2f(box.bottomLeft.X()+box.width, box.bottomLeft.Y());
		glVertex2f(box.bottomLeft.X()+box.width, box.bottomLeft.Y()+box.height);
		glVertex2f(box.bottomLeft.X(), box.bottomLeft.Y()+box.height);
		glVertex2f(box.bottomLeft.X(), box.bottomLeft.Y());
		glEnd();

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glColor3f(box.r, box.g, box.b);
		glBegin(GL_QUADS);
		glVertex2f(box.bottomLeft.X()+box.width, box.bottomLeft.Y());
		glVertex2f(box.bottomLeft.X()+box.width, box.bottomLeft.Y()+box.height);
		glVertex2f(box.bottomLeft.X(), box.bottomLeft.Y()+box.height);
		glVertex2f(box.bottomLeft.X(), box.bottomLeft.Y());
		glEnd();
	}

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glBegin(GL_QUADS);
	if( gAnimationState > 1 ) {
		glColor3f(1, 1, 0);
		Box & box = gBounds;

		glVertex2f(box.bottomLeft.X()+box.width, box.bottomLeft.Y());
		glVertex2f(box.bottomLeft.X()+box.width, box.bottomLeft.Y()+box.height);
		glVertex2f(box.bottomLeft.X(), box.bottomLeft.Y()+box.height);
		glVertex2f(box.bottomLeft.X(), box.bottomLeft.Y());
	}
	glEnd();

	if( gComputeBounds == gAnimationState && gCurrentBox < gBoxCount) {
		Box & box = gBoxes[gCurrentBox];

		glLineWidth(3.0);
		glColor3f(1, 1, 1);
		glBegin(GL_LINES);
		if( gBoxEdge == 0 ) {
			glVertex2f(box.bottomLeft.X(), box.bottomLeft.Y());
			glVertex2f(box.bottomLeft.X(), box.bottomLeft.Y()+box.height);
		} else if( gBoxEdge == 1 ) {
			glVertex2f(box.bottomLeft.X(), box.bottomLeft.Y());
			glVertex2f(box.bottomLeft.X()+box.width, box.bottomLeft.Y());
		} else if( gBoxEdge == 2 ) {
			glVertex2f(box.bottomLeft.X()+box.width, box.bottomLeft.Y());
			glVertex2f(box.bottomLeft.X()+box.width, box.bottomLeft.Y()+box.height);
		} else if( gBoxEdge == 3 ) {
			glVertex2f(box.bottomLeft.X(), box.bottomLeft.Y()+box.height);
			glVertex2f(box.bottomLeft.X()+box.width, box.bottomLeft.Y()+box.height);
		}
		glEnd();
		glLineWidth(1.0);

	}


	glLoadIdentity();

	glEnable(GL_LINE_STIPPLE);
	glLineStipple(2, 0x0f0f);
	glColor3f(1, 1, 1);
	glBegin(GL_LINES); 
	glVertex2i(0, 0);
	glVertex2i(gWindowWidth, 0);
	
	glVertex2i(0, 0);
	glVertex2i(-gWindowWidth, 0);
	
	glVertex2i(0, 0);
	glVertex2i(0, gWindowHeight);
 	
	glVertex2i(0, 0);
	glVertex2i(0, -gWindowHeight);
 	
	glEnd();
	glDisable(GL_LINE_STIPPLE);


	int y = gWindowHeight/2;
	int x = -gWindowWidth/2;
	
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glColor4f(0, 0, 0, .4f);
	glEnable(GL_BLEND);
	glBegin(GL_QUADS);
	glVertex2i(x, y);
	glVertex2i(x, y-3*15);
	glVertex2i(x+220, y-3*15);
	glVertex2i(x+220, y);
	glEnd();
	glDisable(GL_BLEND);

	y -= 15;
	glColor3f(1, 1, 0);
	if( 1 == gAnimationState ) glColor3f(1, 1, 0);
		else glColor3f(1, 1, 1);
		
	glRasterPos2i(x, y);
	{
		std::string text = "1) Rotate Boxes";
		for( int c=0; c<(int)text.size(); ++c ) {
			glutBitmapCharacter(GLUT_BITMAP_9_BY_15, text[c]);
		}
	}
	y -= 15;
	
	if( 2 == gAnimationState ) glColor3f(1, 1, 0);
	else glColor3f(1, 1, 1);
	glRasterPos2i(x, y);
	{
		std::string text = "2) Find x and y extents";
		for( int c=0; c<(int)text.size(); ++c ) {
			glutBitmapCharacter(GLUT_BITMAP_9_BY_15, text[c]);
		}
	}
	y -= 15;
	if( 3 == gAnimationState ) glColor3f(1, 1, 0);
	else glColor3f(1, 1, 1);
	glRasterPos2i(x, y);
	{
		std::string text = "3) Rotate Back";
		for( int c=0; c<(int)text.size(); ++c ) {
			glutBitmapCharacter(GLUT_BITMAP_9_BY_15, text[c]);
		}
	}



	glutSwapBuffers();
}

bool animate = true;
int startTime;
void IdleFunc() {
	static int setCount = 4;

	double speedScale = (setCount-4)/13.0 * 10;

	if( animate ) {
		switch( gAnimationState ) {
			case gRandomize: {
				static bool first = true;

				if( first ) {
					first = false;

					if( ++setCount > 13 ) setCount = 13;
					gBoxCount = max(4, RandomInt(setCount-2, setCount));
					for( int i=0; i<gBoxCount; ++i ) {
						Box & box = gBoxes[i];

						box.bottomLeft.X() = RandomInt(-gWindowWidth*.375, gWindowWidth*.375);
						box.bottomLeft.Y() = RandomInt(-gWindowHeight*.375, gWindowHeight*.375);
						box.width = RandomInt(10, gWindowWidth*.1);
						box.height= RandomInt(10, gWindowHeight*.1);

						box.r = RandomInt(0, 255)/255.0;
						box.g  = RandomInt(0, 255)/255.0;
						box.b = RandomInt(0, 255)/255.0;

					}
				}

				double t = (GetTickCount()-startTime) / 500.0;
				if( t >= 1.0 ) {
					animate = false;
					first = true;
				}

				break;
			}
			case gDeRotate: {
				double t = (GetTickCount()-startTime) / 1000.0;

				if( t >= 1 ) {
					t = 1;
					animate = false;

					gBounds.bottomLeft = gBoxes[0].bottomLeft;
					gBounds.width = 0;
					gBounds.height = 0;
				}
				gAngle = (1-t)*gOriginalAngle + t*0;
				break;
			}
			case gComputeBounds: {
				static bool pause = false;
				static bool first = true;
				if( first ) {
					first = false;
					pause = false;
					gCurrentBox = 0;
					gBoxEdge = 0;					
				}

				Box & box = gBoxes[gCurrentBox];

				if( !pause ) {
					pause = true;

					if( 0 == gBoxEdge ) {
						if( box.bottomLeft.X() < gBounds.bottomLeft.X() ) {
							gBounds.width += abs(gBounds.bottomLeft.X() - box.bottomLeft.X());
							gBounds.bottomLeft.X() = box.bottomLeft.X();
						}

						if( box.bottomLeft.X() > gBounds.bottomLeft.X()+gBounds.width ) {
							gBounds.width += box.bottomLeft.X() - (gBounds.bottomLeft.X()+gBounds.width);
						}
						
					} else if( 1 == gBoxEdge ) {
						if( box.bottomLeft.Y() < gBounds.bottomLeft.Y() ) {
							gBounds.height += gBounds.bottomLeft.Y()-box.bottomLeft.Y();
							gBounds.bottomLeft.Y() = box.bottomLeft.Y();
						}

						if( box.bottomLeft.Y() > gBounds.bottomLeft.Y()+gBounds.height ) {
							gBounds.height += box.bottomLeft.Y() - (gBounds.bottomLeft.Y()+gBounds.height);
						}
					} else if( 2 == gBoxEdge ) {
						if( box.bottomLeft.X()+box.width > gBounds.bottomLeft.X()+gBounds.width ) {
							gBounds.width += (box.bottomLeft.X()+box.width) - (gBounds.bottomLeft.X()+gBounds.width);
						}

						if( box.bottomLeft.X()+box.width < gBounds.bottomLeft.X() ) {
							gBounds.width += gBounds.bottomLeft.X() - (box.bottomLeft.X()+box.width);
							gBounds.bottomLeft.X() = box.bottomLeft.X()+box.width;
						}

					} else if( 3 == gBoxEdge ) {
						if( box.bottomLeft.Y()+box.height > gBounds.bottomLeft.Y()+gBounds.height ) {
							gBounds.height += (box.bottomLeft.Y()+box.height) - (gBounds.bottomLeft.Y()+gBounds.height);
						}

						if( box.bottomLeft.Y()+box.height < gBounds.bottomLeft.Y() ) {
							gBounds.height += gBounds.bottomLeft.Y() - (box.bottomLeft.Y()+box.height);
							gBounds.bottomLeft.Y() = box.bottomLeft.Y()+box.height;
						}

					}
				}

				double t = (GetTickCount()-startTime) / 500.0*speedScale;
				if( t >= 1 ) {
					startTime = GetTickCount();
					pause = false;
					if( 3 == gBoxEdge ) {
						if( ++gCurrentBox >= gBoxCount ) {
							animate = false;
							first = true;
						}
					}
					gBoxEdge = (gBoxEdge+1)%4;

				}
				break;
			}
			case gRotateBack: {
				double t = (GetTickCount()-startTime) / 1000.0;

				if( t >= 1 ) {
					t = 1;
					animate = false;
				}
				gAngle = (t)*gOriginalAngle + t*0;
				break;
			}
			case gDone: {
				double t = (GetTickCount()-startTime) / 500.0;

				if( t >= 1 ) {
					animate = false;
				}
				
			}
		}
		glutPostRedisplay();
	}

	if( !animate ) {
		startTime = GetTickCount();
		animate = true;
		
		if( ++gAnimationState > gDone ) {
			gAnimationState = 0;
		}
	}
}

	

void ReshapeFunc(GLint width, GLint height) {
	gWindowWidth = width;
	gWindowHeight = height;

	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-width/2, width-(width/2), -height/2, height-(height/2), -1, 1);
	glMatrixMode(GL_MODELVIEW);
}

void MouseFunc(int button, int state, int x, int y) {
}

void MotionFunc(int x, int y) {
}

void KeyboardFunc(unsigned char ch, int x, int y) {
	if( 'q' == ch ) {
		exit(0);
	}
}

int main(int argc, char** argv) {

	// GLUT Window Initialization:
	glutInit (&argc, argv);
	glutInitWindowSize (640, 480);
	glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );
	glutCreateWindow ("Rectangle Bounds Detection Visualization");

	gWindowWidth = 640;
	gWindowHeight = 480;

	// Initialize OpenGL graphics state
	InitGraphics();
	InitData();


	// Register callbacks:
	glutDisplayFunc(DisplayFunc);
	glutReshapeFunc(ReshapeFunc);
	glutMouseFunc(MouseFunc);
	glutMotionFunc(MotionFunc);
	glutKeyboardFunc(KeyboardFunc);
	glutIdleFunc(IdleFunc);
	glutMainLoop();

	return 0;
}
 
