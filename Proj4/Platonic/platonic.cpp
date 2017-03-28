/*

Adapted from http://www.cs.umbc.edu/~squire/reference/polyhedra.shtml

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


double *** g_platonicVertexPts; 
int g_kind;

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

void CreateNFF(int nKind, float ra);

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

/*
nKind == 1 ==> Tetrahedron
nKind == 2 ==> Cube
nKind == 3 ==> octahedron
nKind == 4 ==> dodecahedron
nKind == 5 ==> icosahedron
*/
void DrawPlatonic(int nKind)
{
	float variance = 0.5;
	glColor3f(0.9, 0.5, 0);
	if(nKind == 1) // Tetrahedron
	{// 4 poly , 3 vertices per poly
		glBegin(GL_TRIANGLES);
		for(int poly = 0 ; poly < 4 ; poly++)
		{
			for(int vertex = 0 ; vertex < 3 ; vertex++)
			{
				glVertex3f(g_platonicVertexPts[poly][vertex][0],
					g_platonicVertexPts[poly][vertex][1], g_platonicVertexPts[poly][vertex][2]);
				glColor3f(1.0 - variance * vertex, 0.7 - variance * 0.3 * vertex, 0 + variance * vertex *0.2);
			}
		}
		glEnd();
	}
	else if(nKind == 2)	// Cube
	{	// 6 poly, 4 vertices per poly
		glBegin(GL_QUADS);
		for(int poly = 0 ; poly < 6 ; poly++)
			for(int vertex = 0 ; vertex < 4 ; vertex++)
			{
				glVertex3f(g_platonicVertexPts[poly][vertex][0],
					g_platonicVertexPts[poly][vertex][1], g_platonicVertexPts[poly][vertex][2]);
				glColor3f(1.0 - variance * vertex, 0.7 - variance * 0.3 * vertex, 0 + variance * vertex *0.2);
			}
		glEnd();
	}
	else if(nKind == 3) // octahedron
	{	// 8 poly, 3 vertices per poly.
		glBegin(GL_TRIANGLES);
		for(int poly = 0 ; poly < 8 ; poly++)
			for(int vertex = 0 ; vertex < 3 ; vertex++)
			{
				glVertex3f(g_platonicVertexPts[poly][vertex][0],
					g_platonicVertexPts[poly][vertex][1], g_platonicVertexPts[poly][vertex][2]);
				glColor3f(1.0 - variance * vertex, 0.7 - variance * 0.3 * vertex, 0 + variance * vertex *0.2);
			}
		glEnd();
	}
	else if(nKind == 4) // dodecahedron
	{	// 12 poly , 5 vertices per poly.
		for(int poly = 0 ; poly < 12 ; poly++)
		{
			glBegin(GL_POLYGON);
			for(int vertex = 0 ; vertex < 5 ; vertex++)
			{
				glVertex3f(g_platonicVertexPts[poly][vertex][0],
					g_platonicVertexPts[poly][vertex][1], g_platonicVertexPts[poly][vertex][2]);
				glColor3f(1.0 - variance * vertex * 0.2, 0.7 - variance * 0.1 * vertex, 0 + variance * vertex *0.2);
			}
			glEnd();
		}
	}
	else if(nKind == 5) //icosahedron
	{
		for(int poly = 0 ; poly < 20 ; poly++)
		{
			glBegin(GL_TRIANGLES);
			for(int vertex = 0 ; vertex < 3 ; vertex++)
			{
				glVertex3f(g_platonicVertexPts[poly][vertex][0],
					g_platonicVertexPts[poly][vertex][1], g_platonicVertexPts[poly][vertex][2]);
				glColor3f(1.0 - variance * vertex * 0.2, 0.7 - variance * 0.1 * vertex, 0 + variance * vertex *0.2);
			}
			glEnd();
		}
	}
}


void MyDraw()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(1,1,1, 0,0,0, 0,1,0);

	DrawPlatonic(g_kind);

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

double *** GetTetrahedron(float r);
double *** GetCube(float r);
double *** GetOctahedron(float r);
double *** GetDodecahedron(float r);
double *** GetIcosahedron(float r);

int main(int argc, char * argv[])
{
	string arg1, arg2;
	arg1 = argv[1];
	arg2 = argv[2];
	if(argc != 3)
	{
		cout << "Usage> platonic -i #" << endl;
		cout << "1 : tetrahedron, 2: cube , 3: octahedron, 4: dodecahedron, 5: icosahedron" << endl;
		exit(-3);
	}
	if(arg1 != "-i")
	{
		cout << "It must be -i" << endl;
		exit(-2);
	}
	int platonicKind = atoi(arg2.c_str());
	if(platonicKind == 0) exit(-9);
	// At first,  I will create the nff file now.

	g_kind = platonicKind;
	int radius = 1.0f;
	CreateNFF(platonicKind, radius);

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

double *** GetIcosahedron(float r)
{
	double vertices[12][3]; /* 12 vertices with x, y, z coordinates */
	double Pi = 3.141592653589793238462643383279502884197;
	
	double phiaa  = 26.56505; /* phi needed for generation */
	//double r = 1.0; /* any radius in which the polyhedron is inscribed */
	double phia = Pi*phiaa/180.0; /* 2 sets of four points */
	double theb = Pi*36.0/180.0;  /* offset second set 36 degrees */
	double the72 = Pi*72.0/180;   /* step 72 degrees */
	vertices[0][0]=0.0;
	vertices[0][1]=0.0;
	vertices[0][2]=r;
	vertices[11][0]=0.0;
	vertices[11][1]=0.0;
	vertices[11][2]=-r;
	double the = 0.0;
	for(int i=1; i<6; i++)
	{
		vertices[i][0]=r*cos(the)*cos(phia);
		vertices[i][1]=r*sin(the)*cos(phia);
		vertices[i][2]=r*sin(phia);
		the = the+the72;
	}
	the=theb;
	for(int i=6; i<11; i++)
	{
		vertices[i][0]=r*cos(the)*cos(-phia);
		vertices[i][1]=r*sin(the)*cos(-phia);
		vertices[i][2]=r*sin(-phia);
		the = the+the72;
	}

	//map vertices to 20 faces */
	double ***resultVertices = new double**[20];
	for(int i = 0 ; i < 20 ; i++)
	{
		// each polygon has 3 vertices
		resultVertices[i] = new double*[3];
		for(int j = 0 ; j < 3 ; j++)
			// each vertex has 3 coordinates(x,y,z)
			resultVertices[i][j] = new double[3];
	}

	// first polygon => 1-4th vertex => vertex x,y,z (i == 0,1,2)
	for(int i = 0 ; i < 3 ; i++)
	{
		resultVertices[0][0][i] = vertices[0][i];
		resultVertices[0][1][i] = vertices[1][i];
		resultVertices[0][2][i] = vertices[2][i];
	}
	// second polygon
	for(int i = 0 ; i < 3 ; i++)
	{
		resultVertices[1][0][i] = vertices[0][i];
		resultVertices[1][1][i] = vertices[2][i];
		resultVertices[1][2][i] = vertices[3][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{	//3rd polygon
		resultVertices[2][0][i] = vertices[0][i];
		resultVertices[2][1][i] = vertices[3][i];
		resultVertices[2][2][i] = vertices[4][i];
	}

	for(int i = 0 ; i < 3 ; i++)
	{	// 4th polygon
		resultVertices[3][0][i] = vertices[0][i];
		resultVertices[3][1][i] = vertices[4][i];
		resultVertices[3][2][i] = vertices[5][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		// 5th polygon
		resultVertices[4][0][i] = vertices[0][i];
		resultVertices[4][1][i] = vertices[5][i];
		resultVertices[4][2][i] = vertices[1][i];
	}

	for(int i = 0 ; i < 3 ; i++)
	{
		// 6th polygon
		resultVertices[5][0][i] = vertices[11][i];
		resultVertices[5][1][i] = vertices[6][i];
		resultVertices[5][2][i] = vertices[7][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		// 7th polygon
		resultVertices[6][0][i] = vertices[11][i];
		resultVertices[6][1][i] = vertices[7][i];
		resultVertices[6][2][i] = vertices[8][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{	// 8th polygon
		resultVertices[7][0][i] = vertices[11][i];
		resultVertices[7][1][i] = vertices[8][i];
		resultVertices[7][2][i] = vertices[9][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		// 9th polygon
		resultVertices[8][0][i] = vertices[11][i];
		resultVertices[8][1][i] = vertices[9][i];
		resultVertices[8][2][i] = vertices[10][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{	// 10th polygon 11,10,6
		resultVertices[9][0][i] = vertices[11][i];
		resultVertices[9][1][i] = vertices[10][i];
		resultVertices[9][2][i] = vertices[6][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{	// 11th polygon 1,2,6
		resultVertices[10][0][i] = vertices[1][i];
		resultVertices[10][1][i] = vertices[2][i];
		resultVertices[10][2][i] = vertices[6][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{	// 12th polygon (15,16,17,18,19);
		resultVertices[11][0][i] = vertices[2][i];
		resultVertices[11][1][i] = vertices[3][i];
		resultVertices[11][2][i] = vertices[7][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{	// 13th polygon (15,16,17,18,19);
		resultVertices[12][0][i] = vertices[3][i];
		resultVertices[12][1][i] = vertices[4][i];
		resultVertices[12][2][i] = vertices[8][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{	// 14th polygon (15,16,17,18,19);
		resultVertices[13][0][i] = vertices[4][i];
		resultVertices[13][1][i] = vertices[5][i];
		resultVertices[13][2][i] = vertices[9][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{	// 15th polygon (15,16,17,18,19);
		resultVertices[14][0][i] = vertices[5][i];
		resultVertices[14][1][i] = vertices[1][i];
		resultVertices[14][2][i] = vertices[10][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		// 16th polygon (15,16,17,18,19);
		resultVertices[15][0][i] = vertices[6][i];
		resultVertices[15][1][i] = vertices[7][i];
		resultVertices[15][2][i] = vertices[2][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		// 17th polygon (15,16,17,18,19);
		resultVertices[16][0][i] = vertices[7][i];
		resultVertices[16][1][i] = vertices[8][i];
		resultVertices[16][2][i] = vertices[3][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		// 18th polygon (15,16,17,18,19);
		resultVertices[17][0][i] = vertices[8][i];
		resultVertices[17][1][i] = vertices[9][i];
		resultVertices[17][2][i] = vertices[4][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		// 19th polygon (15,16,17,18,19);
		resultVertices[18][0][i] = vertices[9][i];
		resultVertices[18][1][i] = vertices[10][i];
		resultVertices[18][2][i] = vertices[5][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		// 20th polygon (15,16,17,18,19);
		resultVertices[19][0][i] = vertices[10][i];
		resultVertices[19][1][i] = vertices[6][i];
		resultVertices[19][2][i] = vertices[1][i];
	}

	for(int poly = 0 ; poly < 20 ; poly++)
	{
		for(int vertex = 0 ; vertex < 3 ; vertex++)
		{
			resultVertices[poly][vertex][0] /= 5;
			resultVertices[poly][vertex][1] /= 5;
			resultVertices[poly][vertex][2] /= 5;
		}
		glEnd();
	}
	return	resultVertices;
	/*1-polygon(0,1,2);
	2-polygon(0,2,3);
	3-polygon(0,3,4);
	4-polygon(0,4,5);
	5-polygon(0,5,1);
	6-polygon(11,6,7);
	7-polygon(11,7,8);
	8-polygon(11,8,9);
	9-polygon(11,9,10);
	10-polygon(11,10,6);
	11-polygon(1,2,6);
	12-polygon(2,3,7);
	13-polygon(3,4,8);
	14-polygon(4,5,9);
	15-polygon(5,1,10);
	16-polygon(6,7,2);
	17-polygon(7,8,3);
	18-polygon(8,9,4);
	19-polygon(9,10,5);
	20-polygon(10,6,1);*/
}

double *** GetDodecahedron(float r)
{
	double vertices[20][3]; /* 20 vertices with x, y, z coordinate */
	double Pi = 3.141592653589793238462643383279502884197;
	
	double phiaa = 52.62263590; /* the two phi angles needed for generation */
	double phibb = 10.81231754;
	
	//float r = 1.0; /* any radius in which the polyhedron is inscribed */
	float phia = Pi*phiaa/180.0; /* 4 sets of five points each */
	float phib = Pi*phibb/180.0;
	float phic = Pi*(-phibb)/180.0;
	float phid = Pi*(-phiaa)/180.0;
	float the72 = Pi*72.0/180;
	float theb = the72/2.0; /* pairs of layers offset 36 degrees */
	float the = 0.0;
	for(int i=0; i<5; i++)
	{
		vertices[i][0]=r*cos(the)*cos(phia);
		vertices[i][1]=r*sin(the)*cos(phia);
		vertices[i][2]=r*sin(phia);
		the = the+the72;
	}
	the=0.0;
	for(int i=5; i<10; i++)
	{
		vertices[i][0]=r*cos(the)*cos(phib);
		vertices[i][1]=r*sin(the)*cos(phib);
		vertices[i][2]=r*sin(phib);
		the = the+the72;
	}
	the = theb;
	for(int i=10; i<15; i++)
	{
		vertices[i][0]=r*cos(the)*cos(phic);
		vertices[i][1]=r*sin(the)*cos(phic);
		vertices[i][2]=r*sin(phic);
		the = the+the72;
	}
	the=theb;
	for(int i=15; i<20; i++)
	{
		vertices[i][0]=r*cos(the)*cos(phid);
		vertices[i][1]=r*sin(the)*cos(phid);
		vertices[i][2]=r*sin(phid);
		the = the+the72;
	}

	// map vertices to 12 faces(12 poly, 5 vertices)
	double ***resultVertices = new double**[12];
	for(int i = 0 ; i < 12 ; i++)
	{
		// each polygon has 5 vertices
		resultVertices[i] = new double*[5];
		for(int j = 0 ; j < 5 ; j++)
			// each vertex has 3 coordinates(x,y,z)
			resultVertices[i][j] = new double[3];
	}

	// first polygon => 1-4th vertex => vertex x,y,z (i == 0,1,2)
	for(int i = 0 ; i < 3 ; i++)
	{
		resultVertices[0][0][i] = vertices[0][i];
		resultVertices[0][1][i] = vertices[1][i];
		resultVertices[0][2][i] = vertices[2][i];
		resultVertices[0][3][i] = vertices[3][i];
		resultVertices[0][4][i] = vertices[4][i];
	}
	// second polygon
	for(int i = 0 ; i < 3 ; i++)
	{
		resultVertices[1][0][i] = vertices[0][i];
		resultVertices[1][1][i] = vertices[1][i];
		resultVertices[1][2][i] = vertices[6][i];
		resultVertices[1][3][i] = vertices[10][i];
		resultVertices[1][4][i] = vertices[5][i];
	}
	
	//3rd polygon
	for(int i = 0 ; i < 3 ; i++)
	{
		resultVertices[2][0][i] = vertices[1][i];
		resultVertices[2][1][i] = vertices[2][i];
		resultVertices[2][2][i] = vertices[7][i];
		resultVertices[2][3][i] = vertices[11][i];
		resultVertices[2][4][i] = vertices[6][i];
	}

	for(int i = 0 ; i < 3 ; i++)
	{
		// 4th polygon
		resultVertices[3][0][i] = vertices[2][i];
		resultVertices[3][1][i] = vertices[3][i];
		resultVertices[3][2][i] = vertices[8][i];
		resultVertices[3][3][i] = vertices[12][i];
		resultVertices[3][4][i] = vertices[7][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		// 5th polygon
		resultVertices[4][0][i] = vertices[3][i];
		resultVertices[4][1][i] = vertices[4][i];
		resultVertices[4][2][i] = vertices[9][i];
		resultVertices[4][3][i] = vertices[13][i];
		resultVertices[4][4][i] = vertices[8][i];
	}

	for(int i = 0 ; i < 3 ; i++)
	{
		// 6th polygon
		resultVertices[5][0][i] = vertices[4][i];
		resultVertices[5][1][i] = vertices[0][i];
		resultVertices[5][2][i] = vertices[5][i];
		resultVertices[5][3][i] = vertices[14][i];
		resultVertices[5][4][i] = vertices[9][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		// 7th polygon
		resultVertices[6][0][i] = vertices[15][i];
		resultVertices[6][1][i] = vertices[16][i];
		resultVertices[6][2][i] = vertices[11][i];
		resultVertices[6][3][i] = vertices[6][i];
		resultVertices[6][4][i] = vertices[10][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		// 8th polygon
		resultVertices[7][0][i] = vertices[16][i];
		resultVertices[7][1][i] = vertices[17][i];
		resultVertices[7][2][i] = vertices[12][i];
		resultVertices[7][3][i] = vertices[7][i];
		resultVertices[7][4][i] = vertices[11][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		// 9th polygon
		resultVertices[8][0][i] = vertices[17][i];
		resultVertices[8][1][i] = vertices[18][i];
		resultVertices[8][2][i] = vertices[13][i];
		resultVertices[8][3][i] = vertices[8][i];
		resultVertices[8][4][i] = vertices[12][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		// 10th polygon 18,19,14,9,13
		resultVertices[9][0][i] = vertices[18][i];
		resultVertices[9][1][i] = vertices[19][i];
		resultVertices[9][2][i] = vertices[14][i];
		resultVertices[9][3][i] = vertices[9][i];
		resultVertices[9][4][i] = vertices[13][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		// 11th polygon (19,15,10,5,14);
		resultVertices[10][0][i] = vertices[19][i];
		resultVertices[10][1][i] = vertices[15][i];
		resultVertices[10][2][i] = vertices[10][i];
		resultVertices[10][3][i] = vertices[5][i];
		resultVertices[10][4][i] = vertices[14][i];
	}
	for(int i = 0 ; i < 3 ; i++)
	{
		// 12th polygon (15,16,17,18,19);
		resultVertices[11][0][i] = vertices[15][i];
		resultVertices[11][1][i] = vertices[16][i];
		resultVertices[11][2][i] = vertices[17][i];
		resultVertices[11][3][i] = vertices[18][i];
		resultVertices[11][4][i] = vertices[19][i];
	}

	for(int poly = 0 ; poly < 12 ; poly++)
	{
		for(int vertex = 0 ; vertex < 5 ; vertex++)
		{
			resultVertices[poly][vertex][0] /= 5;
			resultVertices[poly][vertex][1] /= 5;
			resultVertices[poly][vertex][2] /= 5;
		}
	}
	return	resultVertices;
	/*1-polygon(0,1,2,3,4);
	
	2-polygon(0,1,6,10,5);
	3-polygon(1,2,7,11,6);
	4-polygon(2,3,8,12,7);
	5-polygon(3,4,9,13,8);
	6-polygon(4,0,5,14,9);
	
	7-polygon(15,16,11,6,10);
	8-polygon(16,17,12,7,11);
	9-polygon(17,18,13,8,12);

	10-polygon(18,19,14,9,13);
	11-polygon(19,15,10,5,14);
	12-polygon(15,16,17,18,19);*/
}

double *** GetOctahedron(float r)
{
	double vertices[6][3]; /* 6 vertices with x, y, z coordinate */
	double Pi = 3.141592653589793238462643383279502884197;
	
	double phiaa  = 0.0; /* the phi needed for generation */
	//float r = 1.0; /* any radius in which the polyhedron is inscribed */
	float phia = Pi*phiaa/180.0; /* 1 set of four points */
	float the90 = Pi*90.0/180;
	vertices[0][0]=0.0;
	vertices[0][1]=0.0;
	vertices[0][2]=r;
	
	vertices[5][0]=0.0;
	vertices[5][1]=0.0;
	vertices[5][2]=-r;
	
	float the = 0.0;
	for(int i=1; i<5; i++)
	{
		vertices[i][0]=r*cos(the)*cos(phia);
		vertices[i][1]=r*sin(the)*cos(phia);
		vertices[i][2]=r*sin(phia);
		the = the+the90;
	}

	// map vertices to 8 faces
	double ***resultVertices = new double**[8];
	for(int i = 0 ; i < 8 ; i++)
	{
		// each polygon has 3 vertices
		resultVertices[i] = new double*[3];
		for(int j = 0 ; j < 3 ; j++)
			// each vertex has 3 coordinates(x,y,z)
			resultVertices[i][j] = new double[3];
	}
	// first polygon => 1-4th vertex => vertex x,y,z (i == 0,1,2)
	for(int i = 0 ; i < 3 ; i++)
	{
		resultVertices[0][0][i] = vertices[0][i];
		resultVertices[0][1][i] = vertices[1][i];
		resultVertices[0][2][i] = vertices[2][i];
	}
	// second polygon
	for(int i = 0 ; i < 3 ; i++)
	{
		resultVertices[1][0][i] = vertices[0][i];
		resultVertices[1][1][i] = vertices[2][i];
		resultVertices[1][2][i] = vertices[3][i];
	}
	//3rd polygon
	for(int i = 0 ; i < 3 ; i++)
	{
		resultVertices[2][0][i] = vertices[0][i];
		resultVertices[2][1][i] = vertices[3][i];
		resultVertices[2][2][i] = vertices[4][i];
	}
	// 4th polygon
	for(int i = 0 ; i < 3 ; i++)
	{
		// 4th polygon
		resultVertices[3][0][i] = vertices[0][i];
		resultVertices[3][1][i] = vertices[4][i];
		resultVertices[3][2][i] = vertices[1][i];
	}
	// 5th polygon
	for(int i = 0 ; i < 3 ; i++)
	{
		// 4th polygon
		resultVertices[4][0][i] = vertices[5][i];
		resultVertices[4][1][i] = vertices[1][i];
		resultVertices[4][2][i] = vertices[2][i];
	}
	// 6th polygon
	for(int i = 0 ; i < 3 ; i++)
	{
		// 4th polygon
		resultVertices[5][0][i] = vertices[5][i];
		resultVertices[5][1][i] = vertices[2][i];
		resultVertices[5][2][i] = vertices[3][i];
	}
	// 7th polygon
	for(int i = 0 ; i < 3 ; i++)
	{
		// 4th polygon
		resultVertices[6][0][i] = vertices[5][i];
		resultVertices[6][1][i] = vertices[3][i];
		resultVertices[6][2][i] = vertices[4][i];
	}
	// 8th polygon
	for(int i = 0 ; i < 3 ; i++)
	{
		// 4th polygon
		resultVertices[7][0][i] = vertices[5][i];
		resultVertices[7][1][i] = vertices[4][i];
		resultVertices[7][2][i] = vertices[1][i];
	}
	for(int poly = 0 ; poly < 8 ; poly++)
	{
		for(int vertex = 0 ; vertex < 3 ; vertex++)
		{
			resultVertices[poly][vertex][0] /= 5;
			resultVertices[poly][vertex][1] /= 5;
			resultVertices[poly][vertex][2] /= 5;
		}
	}
	return	resultVertices;

	/*polygon(0,1,2);
	polygon(0,2,3);
	polygon(0,3,4);
	polygon(0,4,1);
	polygon(5,1,2);
	polygon(5,2,3);
	polygon(5,3,4);
	polygon(5,4,1);*/
}

double *** GetTetrahedron(float r)
{
	double vertices[4][3]; /* 4 vertices with  x, y, z coordinate */
	double Pi = 3.141592653589793238462643383279502884197;
	
	double phiaa  = -19.471220333; /* the phi angle needed for generation */
	
	//float r = 1.0; /* any radius in which the polyhedron is inscribed */
	float phia = Pi*phiaa/180.0; /* 1 set of three points */
	float the120 = Pi*120.0/180.0;
	vertices[0][0] = 0.0;
	vertices[0][1] = 0.0;
	vertices[0][2] = r;
	float the = 0.0;
	for(int i=1; i<4; i++)
	{
	vertices[i][0]=r*cos(the)*cos(phia);
	vertices[i][1]=r*sin(the)*cos(phia);
	vertices[i][2]=r*sin(phia);
	the = the+the120;
	}
	// map vertices to 4 faces
	double ***resultVertices = new double**[4];
	for(int i = 0 ; i < 4 ; i++)
	{
		// each polygon has 3 vertices
		resultVertices[i] = new double*[3];
		for(int j = 0 ; j < 3 ; j++)
			// each vertex has 3 coordinates(x,y,z)
			resultVertices[i][j] = new double[3];
	}
	// first polygon => 1-4th vertex => vertex x,y,z (i == 0,1,2)
	for(int i = 0 ; i < 3 ; i++)
	{
		resultVertices[0][0][i] = vertices[0][i];
		resultVertices[0][1][i] = vertices[1][i];
		resultVertices[0][2][i] = vertices[2][i];
	}
	// second polygon
	for(int i = 0 ; i < 3 ; i++)
	{
		resultVertices[1][0][i] = vertices[0][i];
		resultVertices[1][1][i] = vertices[2][i];
		resultVertices[1][2][i] = vertices[3][i];
	}
	
	//3rd polygon
	for(int i = 0 ; i < 3 ; i++)
	{
		resultVertices[2][0][i] = vertices[0][i];
		resultVertices[2][1][i] = vertices[3][i];
		resultVertices[2][2][i] = vertices[1][i];
	}
	
	for(int i = 0 ; i < 3 ; i++)
	{
		// 4th polygon
		resultVertices[3][0][i] = vertices[1][i];
		resultVertices[3][1][i] = vertices[2][i];
		resultVertices[3][2][i] = vertices[3][i];
	}
	for(int poly = 0 ; poly < 4 ; poly++)
	{
		for(int vertex = 0 ; vertex < 3 ; vertex++)
		{
			resultVertices[poly][vertex][0] /= 5;
			resultVertices[poly][vertex][1] /= 5;
			resultVertices[poly][vertex][2] /= 5;
		}
	}
	return	resultVertices;
	
	/*
	polygon(0,1,2);
	polygon(0,2,3);
	polygon(0,3,1);
	polygon(1,2,3);*/
}

double *** GetCube(float r)
{
	//double vertices[8][3]; /* 8 vertices with x, y, z coordinate */
	double **vertices = new double*[8];
	for(int i = 0 ; i < 8 ; i++)
	{
		vertices[i] = new double[3];
	}
	double Pi = 3.141592653589793238462643383279502884197;
	
	double phiaa = 35.264391; /* the phi needed for generation */
	
	//float r = 1.0; /* any radius in which the polyhedron is inscribed */
	float phia = Pi*phiaa/180.0; /* 2 sets of four points */
	float phib = -phia;
	float the90 = Pi*90.0/180.0;
	float the = 0.0;
	for(int i=0; i<4; i++)
	{
		vertices[i][0]=r*cos(the)*cos(phia);
		vertices[i][1]=r*sin(the)*cos(phia);
		vertices[i][2]=r*sin(phia);
		the = the+the90;
	}
	the=0.0;
	for(int i=4; i<8; i++)
	{
		vertices[i][0]=r*cos(the)*cos(phib);
		vertices[i][1]=r*sin(the)*cos(phib);
		vertices[i][2]=r*sin(phib);
		the = the+the90;
	}
	
	/* map vertices to 6 faces */
	// allocate 6 polygons
	double ***resultVertices = new double**[6];
	for(int i = 0 ; i < 6 ; i++)
	{
		// each polygon has 4 vertices
		resultVertices[i] = new double*[4];
		for(int j = 0 ; j < 4 ; j++)
			// each vertex has 3 coordinates(x,y,z)
			resultVertices[i][j] = new double[3];
	}
	// first polygon => 1-4th vertex => vertex x,y,z (i == 0,1,2)
	for(int i = 0 ; i < 3 ; i++)
	{
		resultVertices[0][0][i] = vertices[0][i];
		resultVertices[0][1][i] = vertices[1][i];
		resultVertices[0][2][i] = vertices[2][i];
		resultVertices[0][3][i] = vertices[3][i];
	}
	// second polygon
	for(int i = 0 ; i < 3 ; i++)
	{
		resultVertices[1][0][i] = vertices[4][i];
		resultVertices[1][1][i] = vertices[5][i];
		resultVertices[1][2][i] = vertices[6][i];
		resultVertices[1][3][i] = vertices[7][i];
	}
	
	//3rd polygon
	for(int i = 0 ; i < 3 ; i++)
	{
		resultVertices[2][0][i] = vertices[0][i];
		resultVertices[2][1][i] = vertices[1][i];
		resultVertices[2][2][i] = vertices[5][i];
		resultVertices[2][3][i] = vertices[4][i];
	}
	
	for(int i = 0 ; i < 3 ; i++)
	{
		// 4th polygon
		resultVertices[3][0][i] = vertices[1][i];
		resultVertices[3][1][i] = vertices[2][i];
		resultVertices[3][2][i] = vertices[6][i];
		resultVertices[3][3][i] = vertices[5][i];
	}
	
	for(int i = 0 ; i < 3 ; i++)
	{
		// 5th polygon
		resultVertices[4][0][i] = vertices[2][i];
		resultVertices[4][1][i] = vertices[3][i];
		resultVertices[4][2][i] = vertices[7][i];
		resultVertices[4][3][i] = vertices[6][i];
	}
	
	for(int i = 0 ; i < 3 ; i++)
	{
		// 6th polygon
		resultVertices[5][0][i] = vertices[3][i];
		resultVertices[5][1][i] = vertices[0][i];
		resultVertices[5][2][i] = vertices[4][i];
		resultVertices[5][3][i] = vertices[7][i];
	}
	for(int poly = 0 ; poly < 6 ; poly++)
	{
		for(int vertex = 0 ; vertex < 4 ; vertex++)
		{
			resultVertices[poly][vertex][0] /= 5;
			resultVertices[poly][vertex][1] /= 5;
			resultVertices[poly][vertex][2] /= 5;
		}
	}
	return	resultVertices;
/*	polygon(0,1,2,3);
	polygon(4,5,6,7);
	polygon(0,1,5,4);
	polygon(1,2,6,5);
	polygon(2,3,7,6);
	polygon(3,0,4,7);*/
}

void CreateNFF(int nKind, float radius)
{
	string kindStr;
	double ***vertexPts;
	switch(nKind)
	{
	case 1 :
		kindStr = "tetrahedron.nff";
		vertexPts = GetTetrahedron(radius);
		break;
	case 2 :
		kindStr = "cube.nff";
		vertexPts = GetCube(radius);
		break;
	case 3 :
		kindStr = "octahedron.nff";
		vertexPts = GetOctahedron(radius);
		break;
	case 4 :
		kindStr = "dodecahedron.nff";
		vertexPts = GetDodecahedron(radius);
		break;
	case 5 :
		kindStr = "icosahedron.nff";
		vertexPts = GetIcosahedron(radius);
		break;
	}
	g_platonicVertexPts = vertexPts;
	ofstream out;
	out.open(kindStr.c_str());
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
	// Tetrahedron Polygons
	if(nKind == 1)
	{
		/*polygon(0,1,2);
		polygon(0,2,3);
		polygon(0,3,1);
		polygon(1,2,3);*/
		for(int poly = 0 ; poly < 4 ; poly ++)
		{
			out << "p 3" << endl;
			for(int vertex = 0 ; vertex < 3 ; vertex++)
			{
				out << vertexPts[poly][vertex][0] << " " << vertexPts[poly][vertex][1] <<
					" " << vertexPts[poly][vertex][2] << endl;
			}
		}
	}
	else if(nKind == 2)	// Cube
	{
		for(int poly = 0 ; poly < 6 ; poly ++)
		{
			out << "p 4" << endl;
			for(int vertex = 0 ; vertex < 4 ; vertex++)
			{
				out << vertexPts[poly][vertex][0] << " " << vertexPts[poly][vertex][1] <<
					" " << vertexPts[poly][vertex][2] << endl;
			}
		}
	}
	else if(nKind == 3)
	{	// octahedron - 8 polygons --> 3 vertices per polygon.
		for(int poly = 0 ; poly < 8 ; poly ++)
		{
			out << "p 3" << endl;
			for(int vertex = 0 ; vertex < 3 ; vertex++)
			{
				out << vertexPts[poly][vertex][0] << " " << vertexPts[poly][vertex][1] <<
					" " << vertexPts[poly][vertex][2] << endl;
			}
		}
	}
	else if(nKind == 4)
	{	//dodecahedron => 12 poly - 5 vertices per poly.
		for(int poly = 0 ; poly < 12 ; poly ++)
		{
			out << "p 5" << endl;
			for(int vertex = 0 ; vertex < 5 ; vertex++)
			{
				out << vertexPts[poly][vertex][0] << " " << vertexPts[poly][vertex][1] <<
					" " << vertexPts[poly][vertex][2] << endl;
			}
		}
	}
	else if(nKind == 5)
	{	// Icos
		for(int poly = 0 ; poly < 20 ; poly ++)
		{
			out << "p 3" << endl;
			for(int vertex = 0 ; vertex < 3 ; vertex++)
			{
				out << vertexPts[poly][vertex][0] << " " << vertexPts[poly][vertex][1] <<
					" " << vertexPts[poly][vertex][2] << endl;
			}
		}
	}
}


