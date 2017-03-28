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
#include <time.h>

int NUM_OF_HOLES = 6;

using namespace std;

class FractalVector{
public:
	float x, y, z, r, theta;

	FractalVector (float _x, float _y, float _r, float _theta) {
		x = _x; //Origin x
		y = _y; //Origin y
		z = 0.0f;
		r = _r; //Length
		theta = _theta; // Angle
	}

	float getEndX() { 
		return x + r*cos(theta/57.3);
	}
	
	float getEndY() {
		return y + r*sin(theta/57.3);
	}
	
	void saveVector() {
		glBegin(GL_QUADS);
			//Draw the current vector
			glColor3f(1.0, 0.3, 0.0);
			glVertex3f(x, y, z);
			glVertex3f(getEndX(), getEndY(), z);
			glColor3f(0.1, 0.4, 0.8);
			glVertex3f(getEndX(), getEndY(), z - 0.1);
			glVertex3f(x, y, z - 0.1);
		glEnd();
		//cout << x << " " << y << " TO " << getEndX() << " " << getEndY() << endl;
	}
};

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
	Point &operator =(const MyVector4 & src)
	{
		x = src.p[0];
		y = src.p[1];
		z = src.p[2];
		return *this;
	}
	operator MyVector4()
	{
		MyVector4 myVec;
		myVec.p[0] =  x;
		myVec.p[1] =  y;
		myVec.p[2] =  z;
		return myVec;
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
				
float holeRadius = 0.15/5;
class Triangle
{
public:
	MyVector4 a, b, c;
	Triangle(MyVector4 _a, MyVector4 _b, MyVector4 _c)
	{
//		_a.p[0] += holeRadius * 2; _a.p[1] += holeRadius * 2/3;
//		_b.p[0] -= holeRadius * 2; _b.p[1] += holeRadius * 2/3;
//		_c.p[1] -= holeRadius * 2;
		a = _a; b = _b; c = _c;
		a.p[0] += holeRadius*1.5; a.p[1] += holeRadius*1.5;
		b.p[0] -= holeRadius*1.5; b.p[1] += holeRadius*1.5;
		c.p[1] -= holeRadius*3.0;
	}

	void DotProduct(float a[3], float b[3], float *result)
	{
		*result = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
	}
	float DotProduct(Point a, Point b)
	{
		float result;
		result = a.x * b.x + a.y * b.y + a.z * b.z;
		return	result;
	}
	void CrossProduct(Point M1, Point M2, Point &resultM)
	{
		resultM.x = M1.y * M2.z - M1.z * M2.y;
		resultM.y = M1.z * M2.x - M1.x * M2.z;
		resultM.z = M1.x * M2.y - M1.y * M2.x;
	}
	// Determination of a point in the triangle
	bool SameSide(Point p1, Point p2, Point a, Point b)
	{
		Point cp1, cp2;
		Point b_a = b - a;
		Point p1_a = p1 - a;
//		cout << " " << p1.x << " " << p1.y << " " << p1.z << endl;
//		cout << " " << a.x << " " << a.y << " " << a.z << endl;
//		cout << " " << p1_a.x << " " << p1_a.y << " " << p1_a.z << endl;
		Point p2_a = p2 - a;
		CrossProduct(b_a, p1_a, cp1);
		CrossProduct(b_a, p2_a, cp2);
		if (DotProduct(cp1, cp2) >= 0)
			return true;
		else
			return false;
	}
	bool PointInTriangle(Point p, Point a, Point b, Point c)
	{
		if(SameSide(p,a, b,c) && SameSide(p,b, a,c) && SameSide(p,c, a,b))
			return true;
		else
			return false;
	}
	bool PointInTriangle(Point p)
	{
		Point pA(a.p[0], a.p[1], a.p[2]);
		Point pB(b.p[0], b.p[1], b.p[2]);
		Point pC(c.p[0], c.p[1], c.p[2]);
		if(SameSide(p, pA, pB, pC) && SameSide(p, pB, pA, pC) && SameSide(p, pC, pA, pB))
			return true;
		else
			return false;
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

MyVector4 cheese[2][3];
MyVector4 cheeseSide[3][4];
MyVector4 hole[12];
//float holeRadius = 0.15;
void MyDraw()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0,0,3, 0,0,0, 0,1,0);
//234 193 82
	glColor3f(234.0/255, 193.0/255, 82.0/255);
	for(int i = 0 ; i < 2 ; i++)
	{
		glBegin(GL_TRIANGLES);
			glVertex3fv(cheese[i][0].p);
			glVertex3fv(cheese[i][1].p);
			glVertex3fv(cheese[i][2].p);
		glEnd();
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		glBegin(GL_QUADS);
			glVertex3fv(cheeseSide[i][0].p);
			glVertex3fv(cheeseSide[i][1].p);
			glVertex3fv(cheeseSide[i][2].p);
			glVertex3fv(cheeseSide[i][3].p);
		glEnd();
	}
//	for(int side = 0 ; side < 2 ; side++)
//	{
		for(int i = 0 ; i < NUM_OF_HOLES ; i++)
		{
			// hole color = 230 162 0
			drawCircle(hole[i].p[0], hole[i].p[1],hole[i].p[2], holeRadius, 230.0/255,162.0/255, 0);
		}
//	}
/*
	for(int i = 0 ; i < 6 ; i+= 3)
	{
		glBegin(GL_TRIANGLES);
			glVertex3fv(g_flower_pole[i].p);
			glVertex3fv(g_flower_pole[i+1].p);
			glVertex3fv(g_flower_pole[i+2].p);
		glEnd();
	}
	glColor3f(1.0, 0.0, 0.0);
	for(int i = 0 ; i < 8 ; i++)
	{
		glBegin(GL_TRIANGLES);
			glVertex3fv(g_flowerbase_left[i].p);
			glVertex3fv(g_flowerbase_right[i].p);
			glVertex3f(0.0f, 0.0f, 0.0f);
		glEnd();
	}
	drawCircle(0, 0, 0, 0.7, 0.5, 0, 0.5);
*/
	glFlush();
}
/*
swisscheese color -> 0.17254902 0.654901961  0.917647059
hole - > 0.164705882 1 0.917647059
MyVector4 cheese[2][3];
MyVector4 cheeseSide[3][4];
MyVector4 hole[2][6];
float holeRadius = 0.2;
*/
void MyKeyboard(unsigned char key, int x, int y)
{
	switch(key)
	{
		case 27 : 
			exit(3);
	}
}

void GetSwisscheese(void);

int main(int argc, char * argv[])
{
	if(argc != 2)
	{
		cout << "USAGE > swisscheese <number_of_holes_In_the_cheese>" << endl;
		exit(-13);
	}
	int nH = atoi(argv[1]);
	NUM_OF_HOLES = nH;
	if(NUM_OF_HOLES > 13)
	{
		cout << "more than 14 holes are crazy!!!!!! " << endl;
		exit(-15);
	}

	GetSwisscheese();
	// At first,  I will create the nff file now.
	CreateNFF("swisscheese.nff");

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

void GetSwisscheese(void)
{
	cheese[0][0] = MyVector4(-1.0/5, 0, 0);	// width = 2, height = 3
	cheese[0][1] = MyVector4( 1.0/5, 0, 0);
	cheese[0][2] = MyVector4( 0, 3.0/5, 0);
	cheese[1][0] = MyVector4(-1.0/5, 0, -1.0/5);	// width = 2, height = 3
	cheese[1][1] = MyVector4( 1.0/5, 0, -1.0/5);
	cheese[1][2] = MyVector4( 0, 3.0/5, -1.0/5);

	cheeseSide[0][0] = MyVector4(0, 3.0/5, 0);
	cheeseSide[0][1] = MyVector4(0, 3.0/5, -1.0/5);
	cheeseSide[0][2] = MyVector4(-1.0/5, 0, -1.0/5);
	cheeseSide[0][3] = MyVector4(-1.0/5, 0, 0);

	cheeseSide[1][0] = MyVector4(0, 3.0/5, 0);
	cheeseSide[1][1] = MyVector4(0, 3.0/5, -1.0/5);
	cheeseSide[1][2] = MyVector4(1.0/5, 0, -1.0/5);
	cheeseSide[1][3] = MyVector4(1.0/5, 0, 0);

	cheeseSide[2][0] = MyVector4(-1.0/5, 0, 0);
	cheeseSide[2][1] = MyVector4(-1.0/5, 0, -1.0/5);
	cheeseSide[2][2] = MyVector4(1.0/5, 0, -1.0/5);
	cheeseSide[2][3] = MyVector4(1.0/5, 0, 0);

	srand(time(NULL));
	int posHoleX;
	int posHoleY;
	float fX;
	float fY;
/*
swisscheese color -> 0.17254902 0.654901961  0.917647059
hole - > 0.164705882 1 0.917647059
*/
	int nHoles = NUM_OF_HOLES;
	float minDist = 10000.0f;
	for(int i = 0 ; i < nHoles ; i++)
	{
		//drawCircle(double x, double y, double z, double radius, float colorR,float colorG,float colorB)
		minDist = 10000.0f;
		posHoleX = rand() % 200;
		posHoleX -= 100;
		posHoleY = rand() % 300;
		fX = float(posHoleX) / 100 / 5.0;
		fY = float(posHoleY) / 100 / 5.0;
		//Check the point is outside the triangle
		Triangle cheese_front(cheese[0][0], cheese[0][1], cheese[0][2]);
		bool inTri = cheese_front.PointInTriangle(Point(fX, fY, 0));
		if(inTri == false)
		{
			i--;
			continue;
		}

		// Detect collision between cheese Holes.
		for(int j = 0 ; j < i ; j++)
		{
			//dist(fX, fY, prev);
			float distance = (fX - hole[j].p[0])*(fX - hole[j].p[0])+(fY-hole[j].p[1])*(fY-hole[j].p[1]);
			distance = sqrt(distance);
			if(minDist > distance)
				minDist = distance;
		}
		if(minDist < holeRadius * 2.3)
		{
//			cout << minDist << endl;
			i--;
			continue;
		}
		hole[i] = MyVector4(fX, fY, 0.002);
		//drawCircle(posHoleX, posHoleY, 0.1, holeRadius, 0.164705882, 1, 0.917647059);
	}
	nHoles = NUM_OF_HOLES;
	minDist = 10000.0f;
	for(int i = nHoles ; i < nHoles * 2 ; i++)
	{
		//drawCircle(double x, double y, double z, double radius, float colorR,float colorG,float colorB)
		minDist = 10000.0f;
		posHoleX = rand() % 200;
		posHoleX -= 100;
		posHoleY = rand() % 300;
		fX = float(posHoleX) / 100 / 5.0;
		fY = float(posHoleY) / 100 / 5.0;
		//Check the point is outside the triangle
		Triangle cheese_front(cheese[0][0], cheese[0][1], cheese[0][2]);
		bool inTri = cheese_front.PointInTriangle(Point(fX, fY, 0));
		if(inTri == false)
		{
			i--;
			continue;
		}

		// Detect collision between cheese Holes.
		for(int j = nHoles ; j < i ; j++)
		{
			//dist(fX, fY, prev);
			float distance = (fX - hole[j].p[0])*(fX - hole[j].p[0])+(fY-hole[j].p[1])*(fY-hole[j].p[1]);
			distance = sqrt(distance);
			if(minDist > distance)
				minDist = distance;
		}
		if(minDist < holeRadius * 2.3)
		{
//			cout << minDist << endl;
			i--;
			continue;
		}
		hole[i] = MyVector4(fX, fY, -1.0/5-0.002);
		//drawCircle(posHoleX, posHoleY, 0.1, holeRadius, 0.164705882, 1, 0.917647059);
	}
}

void CreateNFF(string filename)
{
	ofstream out;
	out.open(filename.c_str());
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
	// cheese polygons(front side)
	for(int poly = 0 ; poly < 2 ; poly++)
	{
		out << "p 3" << endl;
		for(int vertex = 0 ; vertex < 3 ; vertex++)
		{
			out << cheese[poly][vertex].p[0] << " " << cheese[poly][vertex].p[1] << " " << cheese[poly][vertex].p[2] << endl;
		}
	}
	// cheese side polygons
	for(int poly = 0 ; poly < 3 ; poly++)
	{
		out << "p 4" << endl;
		for(int vertex = 0 ; vertex < 4 ; vertex++)
		{
			out << cheeseSide[poly][vertex].p[0] << " " << cheeseSide[poly][vertex].p[1] << " " << cheeseSide[poly][vertex].p[2] << endl;
		}
	}
//	for(int side = 0 ; side < 2 ; side++)
//	{
		for(int iHole = 0 ; iHole < NUM_OF_HOLES * 2 ; iHole++)
		{
			out << "s "; // s centerx y z radius
			out << hole[iHole].p[0] << " " << hole[iHole].p[1] << " " << hole[iHole].p[2] << " " << holeRadius << endl;
		}
//	}

	out.close();

/*
MyVector4 cheese[2][3];
MyVector4 cheeseSide[3][4];
MyVector4 hole[2][6];
float holeRadius = 0.2;
*/
}








