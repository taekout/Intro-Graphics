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
#include <list.h>

#define SCALE_FACTOR 0.005

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
	double saveTranslate[3];
	MyMat4x4()
	{
		SetIdentity();
		saveTranslate[0] = 0.0;
		saveTranslate[1] = 0.0;
		saveTranslate[2] = 0.0;
	}
	void SetIdentity()
	{
		for(int i=0; i<4; i++)
			for(int j=0; j<4; j++)
				p[i][j] = 0.0;
		p[0][0] = 1;	p[1][1] = 1;
		p[2][2] = 1;	p[3][3] = 1;
	}
	void PushTranslate()
	{
		for(int i = 0 ; i < 3 ; i++)
			saveTranslate[i] = p[i][3];
	}
	void PopTranslate()
	{
		for(int i = 0 ; i < 3 ; i++)
			p[i][3] = saveTranslate[i];
	}
	void ClearTranslate()
	{
		for(int i = 0 ; i < 3 ; i++)
			p[i][3] = 0.0f;
	}
	void MultiplyMat(MyMat4x4 m)
	{
		float resultM[4][4];
		int i = 0, j = 0 , k = 0;
		for(i = 0 ; i < 4 ; i++)
		{
			for(j = 0 ; j < 4 ; j++)
			{
				resultM[i][j] = 0.0f;
				for(k = 0 ; k < 4 ; k++)
				{
					resultM[i][j] += m.p[i][k] * p[k][j];
				}
			}
		}
		for(i = 0 ; i < 4 ; i++)
		{
			for(j = 0 ; j < 4 ; j++)
			{
				p[i][j] = resultM[i][j];
			}
		} 
	}
	void RotateX(float theta)
	{
//		SetIdentity();
		MyMat4x4 r;
		float rad = theta * 3.141592/180.0f;
		r.p[1][1] = cos(rad);
		r.p[1][2] = -sin(rad);
		r.p[2][1] = sin(rad);
		r.p[2][2] = cos(rad);
		MultiplyMat(r);
	}
	void RotateZ(float theta)
	{
//		SetIdentity();
		MyMat4x4 r;
		float rad = theta * 3.141592/180.0f;
		r.p[0][0] = cos(rad);
		r.p[0][1] = -sin(rad);
		r.p[1][0] = sin(rad);
		r.p[1][1] = cos(rad);
		MultiplyMat(r);
	}
	void RotateY(float theta)
	{
//		SetIdentity();
		MyMat4x4 r;
		float rad = theta * 3.141592/180.0f;
		r.p[0][0] = cos(rad);
		r.p[0][2] = sin(rad);
		r.p[2][0] = -sin(rad);
		r.p[2][2] = cos(rad);
		MultiplyMat(r);
	}
	void Translate(float x, float y, float z)
	{
//		SetIdentity();
		MyMat4x4 t;
		t.p[0][3] = x;
		t.p[1][3] = y;
		t.p[2][3] = z;
		MultiplyMat(t);
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
		_a.p[0] += holeRadius; _a.p[1] += holeRadius;
		_b.p[0] -= holeRadius; _b.p[1] += holeRadius;
		_c.p[1] -= holeRadius;
		a = _a; b = _b; c = _c;
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

void MyInit()
{
	glClearColor(1.0,1.0,1.0,0.0);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	srand(time(NULL));
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
	gluPerspective(45, 1.0, 0.01, 200);
}

GLfloat fibonacci_series[800] = {};
void CalFibonacci(GLfloat *fibonacci)
{
	for(int i = 2 ; i < 800 ; i++)
	{
		fibonacci[i] = fibonacci[i-1] + fibonacci[i-2];
	}
}

typedef struct QuarterSphere
{
	float v[24 * 24][3];
	QuarterSphere()
	{
		for(int i = 0 ; i < 24 * 24 ; i++)
		{
			for(int j = 0 ; j < 3 ; j++)
			{
				v[i][j] = 0.0f;
			}
		}
	}
} QuarterSphere;

list<QuarterSphere> g_lSphere;
void generateQuarterSphere(int scaley, int scalex, GLfloat r, MyMat4x4 multM, int fibonacci_index = 0)
{
	int i, j;
	float **v;
	v = new float*[scalex * scaley];
	for(int i = 0 ; i < scalex * scaley ; i++)
	{
		v[i] = new float[3];
	}
	//GLfloat v[scalex*scaley][3];
	for (i=0; i<scalex; ++i) {
		for (j=0; j<scaley; ++j) {
			if(fibonacci_index != 0)
			{	// get the previous radius and apply.
				GLfloat prev_r;
				prev_r = fibonacci_series[fibonacci_index - 2];	// previous radius.
				v[i*scaley+j][0]=prev_r*cos(j*2*M_PI/scaley)*cos(i*M_PI/(2*scalex));
			}
			else
			{
				v[i*scaley+j][0]=r*cos(j*2*M_PI/scaley)*cos(i*M_PI/(2*scalex));
			}
			v[i*scaley+j][1]=r*sin(i*M_PI/(2*scalex));
			v[i*scaley+j][2]=r*sin(j*2*M_PI/scaley)*cos(i*M_PI/(2*scalex));

			MyVector4 pt(v[i*scaley+j][0], v[i*scaley+j][1], v[i*scaley+j][2]);
			MyVector4 transPt = multM * pt;
			v[i*scaley+j][0]= transPt.p[0] * SCALE_FACTOR;
			v[i*scaley+j][1]= transPt.p[1] * SCALE_FACTOR;
			v[i*scaley+j][2]= transPt.p[2] * SCALE_FACTOR;
		}
	}

	QuarterSphere qs;
	for(i = 0 ; i < scaley * scalex ; i++)
	{
		for(int j = 0 ; j < 3 ; j++)
		{
			qs.v[i][j] = v[i][j];
		}
	}
	g_lSphere.push_back(qs);

// 	glBegin(GL_LINES);
// 		glColor3f(1, 0, 0);
// 		glVertex3f(0, 0, 0);
// 		glVertex3f(3, 0, 0);
// 		glColor3f(0, 1, 0);
// 		glVertex3f(0, 0, 0);
// 		glVertex3f(0, 3, 0);
// 		glColor3f(0, 0, 1);
// 		glVertex3f(0, 0, 0);
// 		glVertex3f(0, 0, 3);
// 	glEnd();
// 	glBegin(GL_QUADS);
//  		for (i=0; i<scalex-1; ++i) {
//  			for (j=0; j<scaley/2; ++j) {
//  				glColor3f(0, 0.3, 0.5);
//  				glVertex3fv(v[i*scaley+j]); // 1
//  				glColor3f(1, 0.2, 0.1);
//  				glVertex3fv(v[i*scaley+(j+1)%scaley]); // 2
//  				glVertex3fv(v[(i+1)*scaley+(j+1)%scaley]); // 3 
//  				glVertex3fv(v[(i+1)*scaley+j]); // 4
//  			}
//  		}
// 	glEnd();
}

void DrawHalfSpheres(void)
{
	glBegin(GL_LINES);
		glColor3f(1, 0, 0);
		glVertex3f(0, 0, 0);
		glVertex3f(50, 0, 0);
		glColor3f(0, 1, 0);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 50, 0);
		glColor3f(0, 0, 1);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, 50);
	glEnd();
	float fRand1, fRand2, fRand3;
	float red = 0.48627451;
	float green = 0.721568627;
	float blue = 0.964705882;
	int tmp = 0;
	for(list<QuarterSphere>::iterator it = g_lSphere.begin() ; it != g_lSphere.end() ; it++, tmp ++)
	{
		QuarterSphere qs = *it;
		glBegin(GL_QUADS);
		if(tmp != g_lSphere.size() - 1)
		{
			for (int i=0; i<24-1; ++i) {
				for (int j=0; j<24/2; ++j) {
					fRand1 = rand() % 100;
					fRand1 /= 1000.0;
					fRand2 = rand() % 100;
					fRand2 /= 1000.0;
					fRand3 = rand() % 100;
					fRand3 /= 1000.0;
					glColor3f(red + fRand1, green + fRand2, blue + fRand3);
//glColor3f(1.0 ,0 ,0);
					glVertex3fv(qs.v[i*24+j]); // 1
					fRand1 = rand() % 100;
					fRand1 /= 1000.0;
					fRand2 = rand() % 100;
					fRand2 /= 1000.0;
					fRand3 = rand() % 100;
					fRand3 /= 1000.0;
					glColor3f(red + fRand1, green + fRand2, blue + fRand3);
//glColor3f(0.0 ,1 ,0);
					glVertex3fv(qs.v[i*24+(j+1)%24]); // 2
					fRand1 = rand() % 100;
					fRand1 /= 1000.0;
					fRand2 = rand() % 100;
					fRand2 /= 1000.0;
					fRand3 = rand() % 100;
					fRand3 /= 1000.0;
					glColor3f(red + fRand1, green + fRand2, blue + fRand3);
//glColor3f(0.0 ,0 ,1);
					glVertex3fv(qs.v[(i+1)*24+(j+1)%24]); // 3
					fRand1 = rand() % 100;
					fRand1 /= 1000.0;
					fRand2 = rand() % 100;
					fRand2 /= 1000.0;
					fRand3 = rand() % 100;
					fRand3 /= 1000.0;
					glColor3f(red + fRand1, green + fRand2, blue + fRand3);
//glColor3f(1 ,0 ,1);
					glVertex3fv(qs.v[(i+1)*24+j]); // 4
				}
			}
		}
		else
		{
			for (int i=0; i<24-1; ++i) {
				for (int j=0; j<24/2; ++j) {
					fRand1 = rand() % 100;
					fRand1 /= 1000.0;
					fRand2 = rand() % 100;
					fRand2 /= 1000.0;
					fRand3 = rand() % 100;
					fRand3 /= 1000.0;
					glColor3f(red + fRand1, green + fRand2, blue + fRand3);
//glColor3f(1.0 ,0 ,0);
					glVertex3fv(qs.v[i*24+j]); // 1
					fRand1 = rand() % 100;
					fRand1 /= 1000.0;
					fRand2 = rand() % 100;
					fRand2 /= 1000.0;
					fRand3 = rand() % 100;
					fRand3 /= 1000.0;
					glColor3f(red + fRand1, green + fRand2, blue + fRand3);
//glColor3f(0.0 ,1 ,0);
					glVertex3fv(qs.v[i*24+(j+1)%24]); // 2
					fRand1 = rand() % 100;
					fRand1 /= 1000.0;
					fRand2 = rand() % 100;
					fRand2 /= 1000.0;
					fRand3 = rand() % 100;
					fRand3 /= 1000.0;
					glColor3f(red + fRand1, green + fRand2, blue + fRand3);
//glColor3f(0.0 ,0 ,1);
					glVertex3fv(qs.v[(i+1)*24+(j+1)%24]); // 3
					fRand1 = rand() % 100;
					fRand1 /= 1000.0;
					fRand2 = rand() % 100;
					fRand2 /= 1000.0;
					fRand3 = rand() % 100;
					fRand3 /= 1000.0;
					glColor3f(red + fRand1, green + fRand2, blue + fRand3);
//glColor3f(1 ,0 ,1);
					glVertex3fv(qs.v[(i+1)*24+j]); // 4
				}
			}
		}
		glEnd();
	}
}

void GenerateSeashell(int nShells)
{
	fibonacci_series[0] = 1.0f * nShells;
	fibonacci_series[1] = 1.0f * nShells;
	CalFibonacci(fibonacci_series);
	MyMat4x4 mat;
	generateQuarterSphere(24, 24, fibonacci_series[0], mat);
	mat.RotateX(-90);
	generateQuarterSphere(24, 24, fibonacci_series[1], mat);////
	int count4 = 0;
	for(int i = 2 ; i < nShells ; i++)
	{
		// Before Rotate Push Translation Matrix
		mat.PushTranslate();
		mat.ClearTranslate();
		mat.RotateX(-90);
		// After Rotate, Pop Translation Matrix
		mat.PopTranslate();
		cout << "I:"<< i << "   count :" << count4 << endl;
		float fTrans = fibonacci_series[i] - fibonacci_series[i-1];
		switch(count4)
		{
		case 0:
			mat.Translate(0, 0, fTrans);
			break;
		case 1:
			mat.Translate(0, fTrans, 0);
			break;
		case 2:
			mat.Translate(0, 0, -fTrans);
			break;
		case 3:
			mat.Translate(0, -fTrans, 0);
			break;
		}
		count4 ++;
		if(count4 == 4)
			count4 = 0;
		cout << fibonacci_series[i] << endl;
		generateQuarterSphere(24, 24, fibonacci_series[i], mat, i);
	}
}

bool g_ShiftDown;
int g_ButtonDown, g_LastX = -1, g_LastY = -1;
double g_AngleX = 0.0, g_AngleY = 0.0, g_AngleZ = 0.0;
double g_PositionX = 0.0f, g_PositionY = 0.0, g_PositionZ = -4.0;

void mouse(int button, int state, int x, int y) {
    // Record the mouse position when a button is pressed.
    g_LastX = x;
    g_LastY = y;
    g_ButtonDown = button;
    g_ShiftDown = (glutGetModifiers() == GLUT_ACTIVE_SHIFT);
}

void motion(int x, int y) {
    int deltaX = x - g_LastX;
    int deltaY = y - g_LastY;

    switch(g_ButtonDown) {
    case GLUT_LEFT_BUTTON:                  // Rotate in the x-y plane
        g_AngleX -= deltaY;
        g_AngleY -= deltaX;
        break;
    case GLUT_RIGHT_BUTTON:                 // Rotate in the x-z plane
        g_AngleX -= deltaY;
        g_AngleZ -= deltaX;
        break;
    case GLUT_MIDDLE_BUTTON:
        if (g_ShiftDown) {                  // Move in and out
            g_PositionZ -= deltaY * 0.1;
        } else {                            // Move left and right
            g_PositionX -= deltaX * 0.1;
            g_PositionY += deltaY * 0.1;
        }
        break;
    }

    g_LastX = x;
    g_LastY = y;
	glutPostRedisplay();
}

void MyDraw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0,0,0.01, 0,0,0, 0,1,0);

	glTranslated(g_PositionX, g_PositionY, g_PositionZ);
	glRotated(g_AngleX, 1.0, 0.0, 0.0);
	glRotated(g_AngleY, 0.0, 1.0, 0.0);
	glRotated(g_AngleZ, 0.0, 0.0, 1.0);

	DrawHalfSpheres();

	glutSwapBuffers();
}

void MyKeyboard(unsigned char key, int x, int y)
{
	switch(key)
	{
		case 27 : 
			exit(3);
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
	for(list<QuarterSphere>:: iterator it = g_lSphere.begin() ; it != g_lSphere.end() ; it++)
	{
		QuarterSphere qs = *it;
		//qs.v[24 * 24][3];
		for (int i=0; i<24-1; ++i) {
			for (int j=0; j<24/2; ++j) {
				out << "p 4"<< endl;
				out << qs.v[i*24+j][0] << " " << qs.v[i*24+j][1] << " " << qs.v[i*24+j][2] << endl; // 1
				out << qs.v[i*24+(j+1)%24][0] << " " << qs.v[i*24+(j+1)%24][1] << " " << qs.v[i*24+(j+1)%24][2] << endl; // 2
				out << qs.v[(i+1)*24+(j+1)%24][0] << " " << qs.v[(i+1)*24+(j+1)%24][1] << " " << qs.v[(i+1)*24+(j+1)%24][2] << endl; // 3
				out << qs.v[(i+1)*24+j][0] << " " << qs.v[(i+1)*24+j][1] << " " << qs.v[(i+1)*24+j][2] << endl; // 4
			}
		}
	}
}


int main(int argc, char * argv[])
{
	// At first,  I will create the nff file now.
	if(argc != 2)
	{
		cout << "Usage seashell <number> " << endl;
		exit(-3);
	} 
	int seashell_depth = atoi(argv[1]);
	if(seashell_depth < 3)
	{
		cout << "I suggest that the parameter should be bigger than 2 " << endl;
		exit(-5);
	}
	if(seashell_depth > 9)
	{
		cout << "Hey hey.. It's going to be too big! " << endl;
		exit(-6);
	}
	GenerateSeashell(seashell_depth);
	CreateNFF("seashell.nff");

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
	glutInitWindowPosition(100,100);
	glutInitWindowSize(400,400);
	glutCreateWindow("NFF");
	
	MyInit();
	glutReshapeFunc(MyReshape);
	glutDisplayFunc(MyDraw);
	glutKeyboardFunc(MyKeyboard);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutMainLoop();

	return	1;
}
