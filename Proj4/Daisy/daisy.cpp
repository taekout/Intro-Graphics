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


void CreateNFF(string filename, MyVector4 *flowerbase_left,MyVector4 *flowerbase_right,
		MyVector4 *flower_pole);

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


// DAISY DRAW !!!
MyVector4 *g_flowerbase_left;
MyVector4 *g_flowerbase_right;
MyVector4 *g_flower_pole;
void MyDraw()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0,0,-2, 0,0,0, 0,1,0);

	glColor3f(0, 1.0, 0);
	for(int i = 0 ; i < 24 ; i+= 3)
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
	drawCircle(0, 0, 0.01, 0.07, 0.5, 0, 0.050);
	drawCircle(0, 0, -0.01, 0.07, 0.5, 0, 0.050);

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
//	class Point base[8];
//	base[0] = Point(-1, 3, 0);
//	base[1] = Point(1, 3, 0);
	string arg1, arg2;
	arg1 = argv[1];
	arg2 = argv[2];
	if(argc != 3)
	{
		cout << "Usage> daisy -p #" << endl;
		exit(-3);
	}
	if(arg1 != "-p")
	{
		cout << "It must be -p" << endl;
		exit(-2);
	}
	int petalno = atoi(arg2.c_str());
	if(petalno > 8)
	{
		cout << "Number of petals must be less than 9 : <EX> :1 ~ 8" <<endl;
		cout << "The reason is that it just does not look good. not aesthetic." << endl;
		exit(-4);
	}
	float degree = 360.0f / petalno;
	MyVector4 flowerbase_left[petalno];
	MyVector4 flowerbase_right[petalno];
	MyVector4 flowerbase_leftpivot(-0.08, 0.3, 0);
	MyVector4 flowerbase_rightpivot(0.08, 0.3, 0);

	MyVector4 flower_pole[6 * 4];
	flower_pole[0] = MyVector4(-0.03, 0, 0.001);
	flower_pole[1] = MyVector4(0.03, 0, 0.001);
	flower_pole[2] = MyVector4(-0.03, 0.3, 0.001);
	flower_pole[3] = MyVector4(0.03, 0.3, 0.001);
	flower_pole[4] = MyVector4(0.03, 0, 0.001);
	flower_pole[5] = MyVector4(-0.03, 0.3, 0.001);
	MyMat4x4 mat;
	mat.Translate(0, 0, 0.06);
	flower_pole[6] = mat * flower_pole[0];
	flower_pole[7] = mat * flower_pole[1];
	flower_pole[8] = mat * flower_pole[2];
	flower_pole[9] = mat * flower_pole[3];
	flower_pole[10] = mat * flower_pole[4];
	flower_pole[11] = mat * flower_pole[5];
	MyMat4x4 mat2;
	mat2.RotateY(90);\
	MyMat4x4 invTrans;
	invTrans.Translate(0, 0, -0.001);
	MyMat4x4 forwardTrans;
	forwardTrans.Translate(0.03, 0, 0.03);
	flower_pole[12] = forwardTrans * mat2 * invTrans * flower_pole[0];
	flower_pole[13] = forwardTrans * mat2 * invTrans * flower_pole[1];
	flower_pole[14] = forwardTrans * mat2 * invTrans * flower_pole[2];
	flower_pole[15] = forwardTrans * mat2 * invTrans * flower_pole[3];
	flower_pole[16] = forwardTrans * mat2 * invTrans * flower_pole[4];
	flower_pole[17] = forwardTrans * mat2 * invTrans * flower_pole[5];

	MyMat4x4 mat3;
	mat3.RotateY(-90);
	forwardTrans.SetIdentity();
	forwardTrans.Translate(-0.03, 0, 0.03);
	flower_pole[18] = forwardTrans * mat3 * invTrans * flower_pole[0];
	flower_pole[19] = forwardTrans * mat3 * invTrans * flower_pole[1];
	flower_pole[20] = forwardTrans * mat3 * invTrans * flower_pole[2];
	flower_pole[21] = forwardTrans * mat3 * invTrans * flower_pole[3];
	flower_pole[22] = forwardTrans * mat3 * invTrans * flower_pole[4];
	flower_pole[23] = forwardTrans * mat3 * invTrans * flower_pole[5];

	MyMat4x4 rot;
	MyMat4x4 trans;
	trans.Translate(0, -0.4, 0);
	glColor3f(0.0, 1.0, 0.0);
	for(int i = 0 ; i < 6*4 ; i+=3)
	{
		flower_pole[i] = trans * flower_pole[i];
		flower_pole[i+1] = trans * flower_pole[i+1];
		flower_pole[i+2] = trans * flower_pole[i+2];
	}
	g_flower_pole = flower_pole;
	trans.SetIdentity();
	trans.Translate(0, 0.15, 0);
	// Base Coordinates
	for(int i = 0 ; i < petalno ; i++)
	{
		rot.RotateZ(360.0f/petalno * i);
		flowerbase_left[i] =  rot *  flowerbase_leftpivot;
		flowerbase_right[i] = rot * flowerbase_rightpivot;
		rot.SetIdentity();
	}
	g_flowerbase_left = flowerbase_left;
	g_flowerbase_right = flowerbase_right;
	// At first,  I will create the nff file now.
	CreateNFF("daisy.nff", flowerbase_left, flowerbase_right, flower_pole);

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

void CreateNFF(string filename, MyVector4 *flowerbase_left,
		MyVector4 *flowerbase_right, MyVector4 *flower_pole)
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
	// Flower Pole First.
	out << "p 3" << endl;
	for(int i = 0 ; i < 3 ; i++)
	{
		out << flower_pole[i].p[0] << " " << flower_pole[i].p[1] << " " 
				<< flower_pole[i].p[2] << endl;
	}
	out << "p 3" << endl;
	for(int i = 3 ; i < 6 ; i++)
	{
		out << flower_pole[i].p[0] << " " << flower_pole[i].p[1] << " " 
				<< flower_pole[i].p[2] << endl;
	}
	out << "p 3" << endl;
	for(int i = 6 ; i < 9 ; i++)
	{
		out << flower_pole[i].p[0] << " " << flower_pole[i].p[1] << " " 
				<< flower_pole[i].p[2] << endl;
	}
	out << "p 3" << endl;
	for(int i = 9 ; i < 12 ; i++)
	{
		out << flower_pole[i].p[0] << " " << flower_pole[i].p[1] << " " 
				<< flower_pole[i].p[2] << endl;
	}
	out << "p 3" << endl;
	for(int i = 12 ; i < 15 ; i++)
	{
		out << flower_pole[i].p[0] << " " << flower_pole[i].p[1] << " " 
				<< flower_pole[i].p[2] << endl;
	}
	out << "p 3" << endl;
	for(int i = 15 ; i < 18 ; i++)
	{
		out << flower_pole[i].p[0] << " " << flower_pole[i].p[1] << " " 
				<< flower_pole[i].p[2] << endl;
	}
	out << "p 3" << endl;
	for(int i = 18 ; i < 21 ; i++)
	{
		out << flower_pole[i].p[0] << " " << flower_pole[i].p[1] << " " 
				<< flower_pole[i].p[2] << endl;
	}
	out << "p 3" << endl;
	for(int i = 21 ; i < 24 ; i++)
	{
		out << flower_pole[i].p[0] << " " << flower_pole[i].p[1] << " " 
				<< flower_pole[i].p[2] << endl;
	}
	// Center base of Daisy flower first.
	for(int i = 0 ; i < 8 ; i++)
	{
		out << "p 3" << endl;
		out << "0 0 0" << endl; // vertex 1 - always 0,0,0
		out << flowerbase_left[i].p[0] << " " << flowerbase_left[i].p[1] << " " << flowerbase_left[i].p[2] << " " << endl;
		out << flowerbase_right[i].p[0] << " " << flowerbase_right[i].p[1] << " " << flowerbase_right[i].p[2] << " " << endl;
	}
	// Draw Information for the Flower Center Circle
	out << "s 0 0 0.01 0.07" << endl;
	out << "s 0 0 -0.01 0.07" << endl;
}


