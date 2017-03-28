/*

I studied the snowflake math and theory from: http://answers.oreilly.com/topic/1383-create-a-monster-fractal-snowflake-using-processing/

*/

#include <cstdlib>

#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#elif defined(_WIN32)
#include "glut.h"
#else
#include <GL/glut.h>
#endif

#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;


class MyVector4
{
public:
	float p[4];
	MyVector4()
	{
		p[0] = 0; p[1] = 0; p[2] = 0; p[3] = 1.0;
	}
	MyVector4(float x, float y, float z)
	{
		p[0] = x; p[1] = y; p[2] = z; p[3] = 1.0;
	}
};

int g_nPoly = 0;
MyVector4 g_Vertices[200000][4];	// 200000 planes and 4 vertices per polygon
int g_nCylinder = 0;
typedef struct Cylinder
{
	float base_x;
	float base_y;
	float base_z;
	float base_radius;
	float apex_x;
	float apex_y;
	float apex_z;
	float apex_radius;
} Cylinder;

Cylinder g_cylinder[200000][2];

class FractalVector{
public:
	float x, y, z, r, theta;

	FractalVector (float _x, float _y, float _r, float _theta) {
		x = _x; //Origin x
		y = _y; //Origin y
		z = 0.0f;
		r = _r; //Length
		theta = _theta; // Angle
		scaleFactor = 0.5;
	}

	float getEndX() { 
		return x + r*cos(theta/57.3);
	}
	
	float getEndY() {
		return y + r*sin(theta/57.3);
	}

	float scaleFactor;
	void drawVector() {
		float endX = getEndX();
		float endY = getEndY();

		glPushMatrix();
		glTranslatef(x * scaleFactor, y * scaleFactor, z);//1
		glColor3f(0.7, 0.3, 0.4);
		GLUquadric *qobj = gluNewQuadric();
		//gluCylinder(qobj, 0.01, 0.01, 0.3, 24, 24);
		glPopMatrix();
		glPushMatrix();
		glTranslatef(endX * scaleFactor, endY * scaleFactor, z);//2
		//glColor3f(0.7, 0.6, 0.5);
		//gluCylinder(qobj, 0.01, 0.01, 0.3, 24, 24);
		glPopMatrix();

		g_cylinder[g_nCylinder][0].base_x = x * scaleFactor;
		g_cylinder[g_nCylinder][0].base_y = y * scaleFactor;
		g_cylinder[g_nCylinder][0].base_z = z;
		g_cylinder[g_nCylinder][0].base_radius = 0.001;
		g_cylinder[g_nCylinder][0].apex_x = x * scaleFactor;
		g_cylinder[g_nCylinder][0].apex_y = y * scaleFactor;
		g_cylinder[g_nCylinder][0].apex_z = z - 0.3;
		g_cylinder[g_nCylinder][0].apex_radius = 0.001;

		g_cylinder[g_nCylinder][1].base_x = endX * scaleFactor;
		g_cylinder[g_nCylinder][1].base_y = endY * scaleFactor;
		g_cylinder[g_nCylinder][1].base_z = z;
		g_cylinder[g_nCylinder][1].base_radius = 0.001;
		g_cylinder[g_nCylinder][1].apex_x = endX * scaleFactor;
		g_cylinder[g_nCylinder][1].apex_y = endY * scaleFactor;
		g_cylinder[g_nCylinder][1].apex_z = z - 0.3;
		g_cylinder[g_nCylinder][1].apex_radius = 0.001;

		g_nCylinder ++;
		g_nPoly++;

		glBegin(GL_QUADS);
			//Draw the current vector
			glColor3f(0.7, 0.3, 0.4);
			glVertex3f(x * scaleFactor, y * scaleFactor, z);//1
			glVertex3f(endX * scaleFactor, endY * scaleFactor, z);//2
// 			//glColor3f(0.7, 0.6, 0.5);
			glVertex3f(endX * scaleFactor, endY * scaleFactor, z - 0.1);//3
			glVertex3f(x * scaleFactor, y * scaleFactor, z - 0.1);//4
		glEnd();
 		g_Vertices[g_nPoly][0] = MyVector4(x * scaleFactor, y * scaleFactor, z);
 		g_Vertices[g_nPoly][1] = MyVector4(endX * scaleFactor, endY * scaleFactor, z);
 		g_Vertices[g_nPoly][2] = MyVector4(endX * scaleFactor, endY * scaleFactor, z - 0.1);
 		g_Vertices[g_nPoly][3] = MyVector4(x * scaleFactor, y * scaleFactor, z - 0.1);
		//cout << x << " " << y << " TO " << getEndX() << " " << getEndY() << endl;
	}
};

class Point{
public:
	float x;
	float y;
	float z;

public:
	Point(){x = -1; y = -1; z = -1;}
	Point(float x, float y, float z)
	{
		this ->x = x; this ->y = y; this ->z = z;
	}
	Point &operator =(const Point & src)
	{
		x = src.x;
		y = src.y;
		z = src.z;
		return *this;
	}
// ADD * operator!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Point & operator *(float t)
	{
		Point *pt = new Point(this ->x, this ->y, this ->z);
		pt ->x = pt ->x * t;
		pt ->y = pt ->y * t;
		pt ->z = pt ->z * t;
		return *pt;
	}
	Point & operator+=(const Point &src)
	{
		x = x + src.x;
		y = y + src.y;
		z = z + src.z;
		return (*this);
	}
	Point & operator-=(const Point &src)
	{
		x = x - src.x;
		y = y - src.y;
		z = z - src.z;
		return	(*this);
	}
	const Point &operator-(const Point &src) const
	{
		return Point(*this) -= src;
	}
	const Point &operator +(const Point &src) const
	{
		return Point(*this) += src;
	}

	Point GetPoint(void)
	{
		Point ret(x, y, z);
		return	ret;
	}
	void NormalizeVector();
	static float GetLengthVector(Point p);
	static void CrossProduct(const Point p1, const Point p2, Point &p3);
	static float DotProduct(const Point p1, const Point p2)
	{
		return	p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
	}
};

class MyMat4x4
{
public:
	float p[4][4];		//	00 01 02 03
						//  10 11 12 13
						//  20 21 22 23
						//  30 31 32 33
	MyMat4x4()
	{
		SetIdentity();
	}
	void SetIdentity()
	{
		for(int i=0; i<4; i++)
			for(int j=0; j<4; j++)
				p[i][j] = 0.0;
		p[0][0] = 1;	p[1][1] = 1;
		p[2][2] = 1;	p[3][3] = 1;
	}
	void RotateZ(float theta)	
	{
		SetIdentity();
		float rad = theta * 3.141592/180.0f;
		p[0][0] = cos(rad);
		p[0][1] = -sin(rad);
		p[1][0] = sin(rad);
		p[1][1] = cos(rad);
	}
	void Translate(float x, float y, float z)
	{
		SetIdentity();
		p[0][3] = x;
		p[1][3] = y;
		p[2][3] = z;
	}
	MyVector4 operator * (MyVector4 in)
	{
		MyVector4 out;
		for(int j=0; j<4; j++) 
		{
			out.p[j] = 0.0;
			for(int i=0; i<4; i++)
				out.p[j] += p[j][i]*in.p[i]; 
		}
		return out;
	}

	MyMat4x4 operator * (MyMat4x4 in)
	{
		MyMat4x4 out;
		for(int i=0; i<4; i++)
			for(int j=0; j<4; j++)
			{
				out.p[i][j] = 0;
				for(int k=0; k<4; k++)
					out.p[i][j]+=
						p[i][k]*in.p[k][j];
			}
		
		return out;
	}

};				

void CreateNFF(string filename);

void MyInit()
{
	glClearColor(1.0,1.0,1.0,0.0);
}

void drawCircle(double x, double y, double z, double radius, float colorR,float colorG,float colorB)
{
	double y1 = y;
	double x1 = x;
	glColor3f(colorR, colorG, colorB);
	glBegin(GL_TRIANGLES);
	for(int i=0;i<361;i++)
	{
		double angle=(float)(((double)i)/57.29577957795135);   
		double x2=x+(radius*(float)sin((double)angle));
		double y2=y+(radius*(float)cos((double)angle));             
		glVertex3d(x,y, z);
		glVertex3d(x1,y1, z);
		glVertex3d(x2,y2, z);
		y1=y2;
		x1=x2;
	}
	glEnd();
}

void MyReshape(int w, int h)
{
	glViewport(0,0,w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45, 1.0, 0.01, 100);
}

int g_depth = 0;
void fractal(FractalVector *v, int N);
void GenerateSnowflake(void)
{
	g_nCylinder = 0;
	FractalVector *seed = new FractalVector(200.0/400 - 0.5,40.0/400 - 0.5,300.0/400,120);
	fractal(seed,g_depth);
	seed = new FractalVector(350.0/400 - 0.5,300.0/400 - 0.5,300.0/400,-120);
	fractal (seed,g_depth);
	seed = new FractalVector(50.0/400 - 0.5,300.0/400 - 0.5,300.0/400,0);
	fractal (seed,g_depth);
}


void MyDraw()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(1,1,1, 0,0,0, 0,1,0);

	GenerateSnowflake();

	glFlush();
}

void MyKeyboard(unsigned char key, int x, int y)
{
	switch(key)
	{
		case 27 : 
			exit(3);
	}
}

int main(int argc, char * argv[])
{
	string arg1, arg2;
	arg1 = argv[1];
	arg2 = argv[2];
	if(argc != 3)
	{
		cout << "Usage> snowflakes -l #" << endl;
		cout << "# : number of subdivision" << endl;
		exit(-4);
	}
	if(arg1 != "-l")
	{
		cout << "It must be -l" << endl;
		exit(-2);
	}

	g_depth = atoi(arg2.c_str());
	if(g_depth >= 8)
	{
		cout << " depth must be less than 8. (0 - 7) because it generates more than 200000 polygons! Too much!" << endl;
		exit(-8);
	}
	// Generate snowflakes.
	GenerateSnowflake();
	// At first,  I will create the nff file now.
	CreateNFF("snowflake.nff");

	glutInit(&argc, argv);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(400,400);
	glutCreateWindow("NFF");
	
	MyInit();
	glutReshapeFunc(MyReshape);
	glutDisplayFunc(MyDraw);
	glutKeyboardFunc(MyKeyboard);
	glutMainLoop();

	return	1;
}

void fractal(FractalVector *v, int N) {
	if (N == 0) {
		v ->drawVector(); //Draw the current vector
	}
	else{
		FractalVector *t1 = new FractalVector(v ->x,v ->y,v ->r/3.0,v ->theta);
		FractalVector *t2 = new FractalVector(t1 ->getEndX(), t1->getEndY(),
							v ->r/3.0, v ->theta + 60.0);
		FractalVector *t3 = new FractalVector(t2 ->getEndX(), t2->getEndY(),
							v ->r/3.0,v ->theta - 60.0);
		FractalVector *t4 = new FractalVector(t3 ->getEndX(), t3->getEndY(),
							v ->r/3.0,v ->theta);
		fractal(t1,N-1); //Recurse
		fractal(t2,N-1); //Recurse
		fractal(t3,N-1); //Recurse
		fractal(t4,N-1); //Recurse
	}
}

void CreateNFF(string filename)
{
	double ***vertexPts;
//	g_platonicVertexPts = vertexPts;
	ofstream out;
	out.open(filename.c_str());
	cout << filename.c_str() << " will be created." << endl;
	if(!out.is_open())
	{
		cout << "snowflake nff file is not opened." << endl;
		exit(-11);
	}
	out << "b" << " " << "0.078 " << "0.361 " << "0.753" << endl;
	out << "v" << endl;
	out << "from 0 0 0" << endl;
	out << "at 0 0 1" << endl;
	out << "up 0 1 0" << endl;
	out << "angle 45" << endl;
	out << "hither 0.05" << endl;
	out << "resolution 512 512" << endl;
	out << "l 2 28 -5" << endl;
	out << "f 1 0.53 0 1 0 100000 0 0" << endl;
	// Now polygons
	// snowflakes polygons
// 	float base_x;
// 	float base_y;
// 	float base_z;
// 	float base_radius;
// 	float apex_x;
// 	float apex_y;
// 	float apex_z;
// 	float apex_radius;
// c
// %g %g %g %g 
// %g %g %g %g 
	for(int i = 0 ; i < g_nCylinder ; i++)
	{
		out << "c" << endl;
		out << g_cylinder[i][0].base_x << " " << g_cylinder[i][0].base_y << " " << g_cylinder[i][0].base_z << " " << g_cylinder[i][0].base_radius << endl;
		out << g_cylinder[i][0].apex_x << " " << g_cylinder[i][0].apex_y << " " << g_cylinder[i][0].apex_z << " " << g_cylinder[i][0].apex_radius << endl;
		out << "c" << endl;
		out << g_cylinder[i][1].base_x << " " << g_cylinder[i][1].base_y << " " << g_cylinder[i][1].base_z << " " << g_cylinder[i][1].base_radius << endl;
		out << g_cylinder[i][1].apex_x << " " << g_cylinder[i][1].apex_y << " " << g_cylinder[i][1].apex_z << " " << g_cylinder[i][1].apex_radius << endl;
	}

	for(int i = 0 ; i < g_nPoly ; i++)
	{
		out << "p 4" << endl;
		for(int j = 0 ; j < 4 ; j++)
		{
			out << g_Vertices[i][j].p[0] << " " << g_Vertices[i][j].p[1] << " "
				<< g_Vertices[i][j].p[2] << endl;
		}
	}
cout << "Polygon Numbers : " << g_nPoly << endl;

}








