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

#define SCALE_FACTOR 0.1

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
	Point & operator+=(const float src)
	{
		x = x + src;
		y = y + src;
		z = z + src;
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

	float * GetPoint(void)
	{
//		Point ret(x, y, z);
//		return	ret;
		float *p = new float[3];
		p[0] = x;
		p[1] = y;
		p[2] = z;
		return p;
	}
	void NormalizeVector()
	{
		double len = GetLengthVector();
		double reciprocalLen = 1.0 / len;
		x *= reciprocalLen;
		y *= reciprocalLen;
		z *= reciprocalLen;
	}
	double GetLengthVector()
	{
		double len = x * x + y * y + z * z;
		len = sqrt(len);
		return len;
	}

	static void CrossProduct(const Point p1, const Point p2, Point &p3);
	static float DotProduct(const Point p1, const Point p2)
	{
		return	p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
	}
};


class MyMat4x4
{
public:
	float p[4][4];		//  00 01 02 03
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
	MyMat4x4 &operator =(const MyMat4x4 & src)
	{
		for(int i = 0 ; i < 4 ; i++)
		{
			for(int j = 0 ; j < 4 ; j++)
			{
				p[i][j] = src.p[i][j];
			}
		}
		return *this;
	}
};
				
float holeRadius = 0.15/5;
class Triangle
{
public:
	Point a, b, c;
	Triangle(Point _a, Point _b, Point _c)
	{
//		_a.p[0] += holeRadius * 2; _a.p[1] += holeRadius * 2/3;
//		_b.p[0] -= holeRadius * 2; _b.p[1] += holeRadius * 2/3;
//		_c.p[1] -= holeRadius * 2;
		a = _a; b = _b; c = _c;
	}

	Triangle(float _a[3], float _b[3], float _c[3])
	{
		a = Point(_a[0], _a[1], _a[2]);
		b = Point(_b[0], _b[1], _b[2]);
		c = Point(_c[0], _c[1], _c[2]);
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
		Point pA(a.x, a.y, a.z);
		Point pB(b.x, b.y, b.z);
		Point pC(c.x, c.y, c.z);
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
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
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
				prev_r = fibonacci_series[5];	// previous radius.
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
	for(list<QuarterSphere>::iterator it = g_lSphere.begin() ; it != g_lSphere.end() ; it++)
	{
		QuarterSphere qs = *it;
		glBegin(GL_QUADS);
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
		glEnd();
	}
}

void GenerateSeashell(int nShells)
{
	fibonacci_series[0] = 1.0f;
	fibonacci_series[1] = 1.0f;
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

typedef struct SpTriangle
{
	Point p[7];	// [6] == mid point
	float N[3];
	SpTriangle()
	{
		for(int i = 0 ; i < 6 ;i++)
		{
			p[i].x = 0.0f;
			p[i].y = 0.0f;
			p[i].z = 0.0f;

			N[0] = 0.0;
			N[1] = 0.0;
			N[2] = 0.0;
		}
	}
} SpTriangle;

typedef struct MY_SPHERE
{
	int nMaxIndex;
	SpTriangle tri[24][24];
	MY_SPHERE()
	{
		nMaxIndex = 0;
	}
} MY_SPHERE;

MY_SPHERE my_sp;

typedef struct TRI
{
	Point v[3];
} TRI;
list<TRI> g_lTriPorcupine;
list<TRI> g_lTriBody;

void GeneratePorcupine(float r, int n, float spherePortion, float needle_length, MyMat4x4 *mat = NULL, MyMat4x4 *transNeedle= NULL)
{
	my_sp.nMaxIndex = n;
	float dtheta = 180.0f/n /180.0f * 3.141592;
	float dtheta2 = dtheta * 2;
	MyMat4x4 affineMat;
	affineMat.SetIdentity();
	MyMat4x4 tr_needle;
	if(mat != NULL)
	{
		affineMat = *mat;
	}
	if(transNeedle != NULL)
	{
		tr_needle = *transNeedle;
	}

//	glBegin(GL_TRIANGLES);
	for(int i = 0 ; i < n ; i++)
	{
		float theta = dtheta * i;
		float x1 = r * cos(theta);
		float h1 = r * sin(theta);
		float x2 = r * cos(theta + dtheta);
		float h2 = r * sin(theta + dtheta);

		for(int j = 0 ; j < n * spherePortion; j++)
		{
			float theta2 = dtheta2 * j;
			float y1 = cos(theta2);
			float z1 = sin(theta2); // The real Z1 = h1 * z1 , The real Z2 = h
			float y2 = cos(theta2 + dtheta2);
			float z2 = sin(theta2 + dtheta2);

			// calculate a normal vector first
			// N = b x a
			// gain a vector = P2(i+1, j) - P1(i, j)
			float a_x = x2 - x1;
			float a_y = h2 * y1 - h1 * y1;
			float a_z = h2 * z1 - h1 * z1;
			// gain b vector = P3(i, j+1) - P1(i, j)
			float b_x = x1 - x1;
			float b_y = h1 * y2 - h1 * y1;
			float b_z = h1 * z2 - h1 * z1;

			//CrossProduct (b, a);
			float N[3];
			N[0] = b_y * a_z - b_z * a_y;
			N[1] = b_z * a_x - b_x * a_z;
			N[2] = b_x * a_y - b_y * a_x;
			float nLength = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
			if(nLength > 0.000000001)
			{
				N[0] /= nLength;
				N[1] /= nLength;
				N[2] /= nLength;
			}

			// (i , j) == (x1, h1 * y1, h1 * z1);
			// (i+1, j) == (x2, h2 * y1, h2 * z1);
			// (i, j+1) == (x1, h1 * y2, h1 * z2);
			// etc...
			// Triangle points - (i,j+1) (i, j) (i+1, j) ==> 
//			glNormal3f(N[0], N[1], N[2]);
my_sp.tri[i][j].N[0] = N[0];
my_sp.tri[i][j].N[1] = N[1];
my_sp.tri[i][j].N[2] = N[2];
//			glColor3f(1, 0, 0);
			Point pt1(x1, h1 * y2, h1 * z2);
			Point pt2(x2, h2 * y1, h2 * z1);
			Point pt3(x1, h1 * y1, h1 * z1);
//			glVertex3f(pt1.x, pt1.y, pt1.z);
//			glVertex3f(pt2.x, pt2.y, pt2.z);
//			glVertex3f(pt3.x, pt3.y, pt3.z);
pt1 = affineMat * pt1;
pt2 = affineMat * pt2;
pt3 = affineMat * pt3;
my_sp.tri[i][j].p[0] = pt1;
my_sp.tri[i][j].p[1] = pt2;
my_sp.tri[i][j].p[2] = pt3;
			// (i+1,j) (i+1, j+1), (i, j+1)
			//glColor3f(0, 0, 1);
			pt1 = Point(x2, h2 * y1, h2 * z1);
			pt2 = Point(x1, h1 * y2, h1 * z2);
			pt3 = Point(x2, h2 * y2, h2 * z2);
//			glVertex3f(pt1.x, pt1.y, pt1.z);
//			glVertex3f(pt2.x, pt2.y, pt2.z);
//			glVertex3f(pt3.x, pt3.y, pt3.z);
pt1 = affineMat * pt1;
pt2 = affineMat * pt2;
pt3 = affineMat * pt3;
my_sp.tri[i][j].p[3] = pt1;
my_sp.tri[i][j].p[4] = pt2;
my_sp.tri[i][j].p[5] = pt3;
Point mid_pt = my_sp.tri[i][j].p[0] + my_sp.tri[i][j].p[2] + my_sp.tri[i][j].p[1] + my_sp.tri[i][j].p[5];
mid_pt = mid_pt * 0.25f * needle_length;
mid_pt = tr_needle * mid_pt;
my_sp.tri[i][j].p[6] = mid_pt;

// Save the vertices for each triangles.
TRI needles;
needles.v[0] = my_sp.tri[i][j].p[6];
needles.v[1] = my_sp.tri[i][j].p[1];
needles.v[2] = my_sp.tri[i][j].p[5];
g_lTriPorcupine.push_back(needles);
needles.v[0] = my_sp.tri[i][j].p[6];
needles.v[1] = my_sp.tri[i][j].p[2];
needles.v[2] = my_sp.tri[i][j].p[1];
g_lTriPorcupine.push_back(needles);
needles.v[0] = my_sp.tri[i][j].p[0];
needles.v[1] = my_sp.tri[i][j].p[2];
needles.v[2] = my_sp.tri[i][j].p[6];
g_lTriPorcupine.push_back(needles);
needles.v[0] = my_sp.tri[i][j].p[5];
needles.v[1] = my_sp.tri[i][j].p[0];
needles.v[2] = my_sp.tri[i][j].p[6];
g_lTriPorcupine.push_back(needles);
		}
	}
//	glEnd();
}

void GeneratePorcupineHead(float r, int n, float needle_length, float trans_needles, MyMat4x4 *mat)
{
	my_sp.nMaxIndex = n;
	float dtheta = 180.0f/n /180.0f * 3.141592;
	float dtheta2 = dtheta * 2;
	MyMat4x4 affineMat;
	affineMat.SetIdentity();
	if(mat != NULL)
	{
		affineMat = *mat;
	}

	for(int i = 0 ; i < n ; i++)
	{
		float theta = dtheta * i;
		float x1 = r * cos(theta);
		float h1 = r * sin(theta);
		float x2 = r * cos(theta + dtheta);
		float h2 = r * sin(theta + dtheta);

		for(int j = 0 ; j < n ; j++)
		{
			float theta2 = dtheta2 * j;
			float y1 = cos(theta2);
			float z1 = sin(theta2); // The real Z1 = h1 * z1 , The real Z2 = h
			float y2 = cos(theta2 + dtheta2);
			float z2 = sin(theta2 + dtheta2);

			// calculate a normal vector first
			// N = b x a
			// gain a vector = P2(i+1, j) - P1(i, j)
			float a_x = x2 - x1;
			float a_y = h2 * y1 - h1 * y1;
			float a_z = h2 * z1 - h1 * z1;
			// gain b vector = P3(i, j+1) - P1(i, j)
			float b_x = x1 - x1;
			float b_y = h1 * y2 - h1 * y1;
			float b_z = h1 * z2 - h1 * z1;

			//CrossProduct (b, a);
			float N[3];
			N[0] = b_y * a_z - b_z * a_y;
			N[1] = b_z * a_x - b_x * a_z;
			N[2] = b_x * a_y - b_y * a_x;
			float nLength = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
			if(nLength > 0.000000001)
			{
				N[0] /= nLength;
				N[1] /= nLength;
				N[2] /= nLength;
			}

			// (i , j) == (x1, h1 * y1, h1 * z1);
			// (i+1, j) == (x2, h2 * y1, h2 * z1);
			// (i, j+1) == (x1, h1 * y2, h1 * z2);
			// etc...
			// Triangle points - (i,j+1) (i, j) (i+1, j) ==> 
//			glNormal3f(N[0], N[1], N[2]);
my_sp.tri[i][j].N[0] = N[0];
my_sp.tri[i][j].N[1] = N[1];
my_sp.tri[i][j].N[2] = N[2];
//			glColor3f(1, 0, 0);
			Point pt1(x1, h1 * y2, h1 * z2);
			Point pt2(x2, h2 * y1, h2 * z1);
			Point pt3(x1, h1 * y1, h1 * z1);
//			glVertex3f(pt1.x, pt1.y, pt1.z);
//			glVertex3f(pt2.x, pt2.y, pt2.z);
//			glVertex3f(pt3.x, pt3.y, pt3.z);
pt1 = affineMat * pt1;
pt2 = affineMat * pt2;
pt3 = affineMat * pt3;
my_sp.tri[i][j].p[0] = pt1;
my_sp.tri[i][j].p[1] = pt2;
my_sp.tri[i][j].p[2] = pt3;
			// (i+1,j) (i+1, j+1), (i, j+1)
			//glColor3f(0, 0, 1);
			pt1 = Point(x2, h2 * y1, h2 * z1);
			pt2 = Point(x1, h1 * y2, h1 * z2);
			pt3 = Point(x2, h2 * y2, h2 * z2);
//			glVertex3f(pt1.x, pt1.y, pt1.z);
//			glVertex3f(pt2.x, pt2.y, pt2.z);
//			glVertex3f(pt3.x, pt3.y, pt3.z);
pt1 = affineMat * pt1;
pt2 = affineMat * pt2;
pt3 = affineMat * pt3;
my_sp.tri[i][j].p[3] = pt1;
my_sp.tri[i][j].p[4] = pt2;
my_sp.tri[i][j].p[5] = pt3;
Point mid_pt = my_sp.tri[i][j].p[0] + my_sp.tri[i][j].p[2] + my_sp.tri[i][j].p[1] + my_sp.tri[i][j].p[5];
mid_pt = mid_pt * 0.25f;	// Get the mean value.
MyMat4x4 inverseTrans;
inverseTrans.Translate(-affineMat.p[0][3], -affineMat.p[1][3], -affineMat.p[2][3]);
mid_pt = inverseTrans * mid_pt;
mid_pt = mid_pt * 2.5;
MyMat4x4 only_trans_back;
only_trans_back.Translate(affineMat.p[0][3], affineMat.p[1][3], affineMat.p[2][3]);
mid_pt = only_trans_back * mid_pt;
MyMat4x4 transNeedleMat;
transNeedleMat.Translate(trans_needles * 1.3, trans_needles , trans_needles);
transNeedleMat.RotateX(80);
mid_pt = transNeedleMat * mid_pt;
my_sp.tri[i][j].p[6] = mid_pt;

// Save the vertices for each triangles.
TRI needles;
needles.v[0] = my_sp.tri[i][j].p[5];
needles.v[1] = my_sp.tri[i][j].p[1];
needles.v[2] = my_sp.tri[i][j].p[6];
g_lTriPorcupine.push_back(needles);
needles.v[0] = my_sp.tri[i][j].p[1];
needles.v[1] = my_sp.tri[i][j].p[2];
needles.v[2] = my_sp.tri[i][j].p[6];
g_lTriPorcupine.push_back(needles);
needles.v[0] = my_sp.tri[i][j].p[6];
needles.v[1] = my_sp.tri[i][j].p[2];
needles.v[2] = my_sp.tri[i][j].p[0];
g_lTriPorcupine.push_back(needles);
needles.v[0] = my_sp.tri[i][j].p[6];
needles.v[1] = my_sp.tri[i][j].p[0];
needles.v[2] = my_sp.tri[i][j].p[5];
g_lTriPorcupine.push_back(needles);
		}
	}
}


void GenerateMySphere(float r, int n)
{
	MyMat4x4 mat;
	mat.Translate(0.0, 0., 0.0);
	GeneratePorcupine(r, n, 1.0, 1.0, &mat);
}


float * ComputeNormal(Triangle tri)
{
	float *N = new float[3];
	Point v1 = tri.b - tri.a;
	Point v2 = tri.c - tri.a;
	Point crossRes;
	tri.CrossProduct(v1, v2, crossRes);
	crossRes.NormalizeVector();
	N[0] = crossRes.x;
	N[1] = crossRes.y;
	N[2] = crossRes.z;

	return	N;
}

void DrawPorcupine(float needle_length)
{
	MyMat4x4 mat;
	mat.RotateX(-125);
	glBegin(GL_TRIANGLES);
		for(list<TRI>::iterator it = g_lTriPorcupine.begin() ; it != g_lTriPorcupine.end() ; it++)
		{
			TRI t = *it;
			float *N;
			Triangle tri_class(t.v[0], t.v[1], t.v[2]);
			N = ComputeNormal(tri_class);
			glNormal3fv(N);
			for(int i = 0 ; i < 3 ; i++)
			{
				MyVector4 v = (MyVector4)t.v[i];
				MyVector4 p = mat * v;
				glVertex3f(p.p[0], p.p[1], p.p[2]);
			}
		}
		for(list<TRI>::iterator it = g_lTriBody.begin() ; it != g_lTriBody.end() ; it++)
		{
			TRI t = *it;
			float *N;
			Triangle tri_class(t.v[0], t.v[1], t.v[2]);
//			N = ComputeNormal(tri_class);
//			glNormal3fv(N);
			for(int i = 0 ; i < 3 ; i++)
			{
				MyVector4 v = (MyVector4)t.v[i];
				MyVector4 p = mat * v;
				glVertex3f(p.p[0], p.p[1], p.p[2]);
			}
		}
	glEnd();
/*	for(int i = 0 ; i < my_sp.nMaxIndex ; i++)
	{
		for(int j = 0 ; j < my_sp.nMaxIndex ; j++)
		{
			float * fp[7];
			fp[0] = my_sp.tri[i][j].p[0].GetPoint();
			fp[1] = my_sp.tri[i][j].p[1].GetPoint();
			fp[2] = my_sp.tri[i][j].p[2].GetPoint();
			fp[3] = my_sp.tri[i][j].p[3].GetPoint();
			fp[4] = my_sp.tri[i][j].p[4].GetPoint();
			fp[5] = my_sp.tri[i][j].p[5].GetPoint();
			Point mid_pt;
			mid_pt = my_sp.tri[i][j].p[0] + my_sp.tri[i][j].p[2] + my_sp.tri[i][j].p[1] + my_sp.tri[i][j].p[5];
			mid_pt = mid_pt * 0.25f * needle_length;
			fp[6] = mid_pt.GetPoint();
// 			glBegin(GL_QUADS);
// glColor3f(0,1,0);
// 				glNormal3fv(my_sp.tri[i][j].N);
// 				glVertex3fv(fp[5]);
// 				glVertex3fv(fp[1]);
// 				glVertex3fv(fp[2]);
// 				glVertex3fv(fp[0]);
// 			glEnd();
			glBegin(GL_TRIANGLES);
				Triangle tri(fp[6], fp[1], fp[5]);
				float *N;
				N = ComputeNormal(tri);
				glNormal3fv(N);
				glVertex3fv(fp[6]);
				glVertex3fv(fp[1]);
				glVertex3fv(fp[5]);

				tri = Triangle(fp[6], fp[2], fp[1]);
				N = ComputeNormal(tri);
				glNormal3fv(N);
				glVertex3fv(fp[6]);
				glVertex3fv(fp[2]);
				glVertex3fv(fp[1]);

				tri = Triangle(fp[0], fp[2], fp[6]);
				N = ComputeNormal(tri);
				glNormal3fv(N);
				glVertex3fv(fp[0]);
				glVertex3fv(fp[2]);
				glVertex3fv(fp[6]);

				tri = Triangle(fp[5], fp[0], fp[6]);
				N = ComputeNormal(tri);
				glNormal3fv(N);
				glVertex3fv(fp[5]);
				glVertex3fv(fp[0]);
				glVertex3fv(fp[6]);
			glEnd();
			for(int k = 0 ; k < 7;k++)
			{
				delete [] fp[k];
			}
		}
	}*/
}

void MyDraw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0,0,2, 0,0,0, 0,1,0);

	glTranslated(g_PositionX, g_PositionY, g_PositionZ);
	glRotated(g_AngleX, 1.0, 0.0, 0.0);
	glRotated(g_AngleY, 0.0, 1.0, 0.0);
	glRotated(g_AngleZ, 0.0, 0.0, 1.0);

	glBegin(GL_LINES);
		glColor3f(1,0,0);
		glVertex3f(0, 0, 0);
		glVertex3f(100, 0, 0);
		glColor3f(0,1,0);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 100, 0);
		glColor3f(0,0,1);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, 100);
	glEnd();

	glColor3f(0, 0.5, 0.5);
	DrawPorcupine(1.8f);

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

	MyMat4x4 mat;
	mat.RotateX(-125);
	for(list<TRI>::iterator it = g_lTriPorcupine.begin() ; it != g_lTriPorcupine.end() ; it++)
	{
		TRI t = *it;
		float *N;
		Triangle tri_class(t.v[0], t.v[1], t.v[2]);
		N = ComputeNormal(tri_class);
		glNormal3fv(N);
		out << "p 3" << endl;
		for(int i = 0 ; i < 3 ; i++)
		{
			Point p = t.v[i];
			p = mat * p;
			out << p.x << " " << p.y << " " << p.z << endl;
			//glVertex3fv(t.v[i].GetPoint());
		}
	}
	for(list<TRI>::iterator it = g_lTriBody.begin() ; it != g_lTriBody.end() ; it++)
	{
		TRI t = *it;
		float *N;
		Triangle tri_class(t.v[0], t.v[1], t.v[2]);
		N = ComputeNormal(tri_class);
		glNormal3fv(N);
		out << "p 3" << endl;
		for(int i = 0 ; i < 3 ; i++)
		{
			Point p = t.v[i];
			p = mat * p;
			out << p.x << " " << p.y << " " << p.z << endl;
			//glVertex3fv(t.v[i].GetPoint());
		}
	}
}

int main(int argc, char * argv[])
{
	// At first,  I will create the nff file now.
	MyMat4x4 transNeedle;
	transNeedle.Translate(0.6, 0.0, 0.0);
	GeneratePorcupine(1, 24, 0.7, 2.3, NULL, &transNeedle);
	GenerateMySphere(0.9, 24);
	MyMat4x4 mat;
	//mat.RotateZ(-90);
	mat.Translate(-1.5, +0.5, 0);
	GeneratePorcupineHead(0.5, 24, 0, 0.5, &mat);
	CreateNFF("porcupine.nff");

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
