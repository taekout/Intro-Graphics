/*
 * Simple GL example
 */

#include "view.h"
#include "motion.h"

// Apple's annoying non-standard GL include location
#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <iostream>
#include <stdio.h>
#include <list>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using namespace std;


typedef struct POINT3
{
	float x;
	float y;
	float z;
} POINT3;

typedef struct VECTOR3
{
	float x;
	float y;
	float z;
} VECTOR3;

typedef struct VIEW
{
	POINT3	camera;
        POINT3	lookAt;
        VECTOR3	cameraUp;
        float	angle;
        float	dNear;		// In this case, the view plane(projection plane) lies at the near Plane
        struct	RESOLUTION
	{
		int x;
		int y;
	} resolution;
} VIEW;

//typedef Point3 Color3;
typedef struct COLOR3
{
	float R;
	float G;
	float B;
} COLOR3;

typedef struct SHADING
{
	COLOR3	color;
	float	Kd;
	float	Ks;
	float	PhongPow;
	float	Transmission;
	float	indexRefraction;
} SHADING;

typedef struct SPHERE
{
	float x;
	float y;
	float z;
	float radius;
} SPHERE;

typedef struct TETRA
{
	float a[3];
	float b[3];
	float c[3];
	float N[3];
	float D[3];
} TETRA;

typedef struct TRIANGLE_STRUCTURE
{
	float a[3];
	float b[3];
	float c[3];
	float d[3];
	float N[3];
	float D[3];
} TRI_STR;

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
	/*Point &operator =(const Point & src)
	{
		x = src.x;
		y = src.y;
		z = src.z;
		return *this;
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
		return Point(*this) -= src;
	}*/

	Point GetPoint(void)
	{
		Point ret(x, y, z);
		return	ret;
	}
	void NormalizeVector();
	static float GetLengthVector(Point p);
	static void CrossProduct(const Point p1, const Point p2, Point &p3);
};

class Parser
{
public:
	ifstream	m_file;

	int m_nTetra;
	list<TRI_STR>	m_lTri;
	Point		m_light;

	VIEW		m_view;
	COLOR3		m_background;
	SHADING		m_shading;

	Point		m_cameraVector;
	void Read(string filename);

	// Member Functions
	void DoViewPoint();
	void DoLightSource();
	void DoBackground();
	void DoFill();
	void DoCone();
	void DoSphere();
	void DoPoly();
	void DoComment();
};


class Viewing{
public:
/*typedef struct VIEW
{
	POINT3	camera;
        POINT3	lookAt;
        VECTOR3	cameraUp;
        float	angle;
        float	dNear;		// In this case, the view plane(projection plane) lies at the near Plane
        struct	RESOLUTION
	{
		int x;
		int y;
	} resolution;
} VIEW;*/
	Viewing()
	{
	}
	~Viewing()
	{
	}

	// Member Functions
	float *ComputeViewMatrix(Point **cameraUp, Point **cameraRight, Point **cameraViewDir);
	void SetView(VIEW view)
	{
		cameraPos = new Point(view.camera.x, view.camera.y, view.camera.z);
		centerPos = new Point(view.lookAt.x, view.lookAt.y, view.lookAt.z);
		UP = new Point(view.cameraUp.x, view.cameraUp.y, view.cameraUp.z);
		this ->angle = view.angle;
		this ->dNear = dNear;
		x_resolution = view.resolution.x;
		y_resolution = view.resolution.y;
		// Compute Camera Vectors
		viewMatrix = ComputeViewMatrix(&m_cameraUp, &m_cameraRight, &m_cameraOppositeViewDir);
	}

// Member Variables
	Point *centerPos;
	Point *cameraPos;
	Point *UP;	// This is just world coordinate uP vector. We need to compute camera vectors separately.
	float angle;
	float dNear;
	float *viewMatrix;
	Point *m_cameraUp, *m_cameraRight, *m_cameraOppositeViewDir;
	int x_resolution, y_resolution;

};

Parser g_parser;
Viewing g_view;
string nff_filename;


void draw(void);

void key(unsigned char k, int x, int y);

/* initialize GLUT - windows and interaction */
void initGLUT(int *argcp, char *argv[])
{
  /* ask for a window at 0,0 with dimensions winWidth, winHeight */
  /* need color, depth (for 3D drawing) and double buffer (smooth display) */
  glutInit(argcp, argv);
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(winWidth, winHeight);
  glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
  glutCreateWindow("Modeled Scene: 435/634 Assignment 2");

  /* set callback functions to be called by GLUT for drawing, window
     resize, keypress, mouse button press, and mouse movement */
  glutDisplayFunc(draw);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(key);
  glutMouseFunc(mousePress);
  glutMotionFunc(mouseDrag);
  g_parser.Read(nff_filename);
}

/* initialize OpenGL - rendering state */
void initGL()
{
  float lightdir[4] = {1,1,2,0};	/* homogeneous light position: directional if w=0 */
  float white[4] = {1,1,1,1}; 		/* color for light: glLightfv needs 4 components!*/
  float dim[4] = {.2,.2,.2,1};
 
  /* enable some GL features */
  glEnable(GL_DEPTH_TEST);		/* tell OpenGL to handle overlapping surfaces */
  glEnable(GL_COLOR_MATERIAL);		/* map glColor to surface colors used by lighting */
  glEnable(GL_NORMALIZE);		/* tell GL to normalize normals so we don't have to */

  /* set up one light for both directional and ambient */
  glLightfv(GL_LIGHT0, GL_AMBIENT, dim);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
  glLightfv(GL_LIGHT0, GL_POSITION, lightdir);
  glEnable(GL_LIGHT0);			/* turn on this light */
  glEnable(GL_LIGHTING);		/* turn on use of lighting in general */
}

int main(int argc, char *argv[])
{ 
  /* set up GLUT and OpenGL */
  initGLUT(&argc, argv);
  nff_filename = "level.nff";
  initGL();
  
  /* let glut take over, it goes into a loop, checking for input and
     calling the input callbacks, then seeing if we need to draw and
     calling the draw callback, ad infinitum */
  glutMainLoop();

  return 0;             /* ANSI C requires main to return int. */
}





void Point::NormalizeVector()
{
	float length = Point::GetLengthVector(*this);
	x /= length;
	y /= length;
	z /= length;
}

float Point::GetLengthVector(Point p)
{
	float length = p.x * p.x + p.y * p.y + p.z * p.z;
	length = sqrt(length);
	return	length;
}

float *Viewing::ComputeViewMatrix(Point **cameraUp, Point **cameraRight, Point **cameraOppositeViewDir)
{
	Point f;
	f.x = g_view.centerPos ->x - g_view.cameraPos ->x; // f == view vector
	f.y = g_view.centerPos ->y - g_view.cameraPos ->y;
	f.z = g_view.centerPos ->z - g_view.cameraPos ->z;
	f.NormalizeVector();
	Point up;
	up.x = g_view.UP ->x;
	up.y = g_view.UP ->y;
	up.z = g_view.UP ->z;
	up.NormalizeVector();

	Point s;
	Point::CrossProduct(f, up, s);	// s == camera right vector
	s.NormalizeVector();
	Point u;
	Point::CrossProduct(s, f, u);	// u == camera up vector
	u.NormalizeVector();

	*cameraUp = new Point(u.x, u.y, u.z);
	*cameraRight = new Point(s.x, s.y, s.z);
	*cameraOppositeViewDir = new Point(-f.x, -f.y, -f.z);
 
	// Check the order !!! -> Vertical or Horizontal?
	float * M = new float[16];
	M[0] = s.x; M[1] = s.y; M[2] = s.z; M[3] = 0;
	M[4] = u.x; M[5] = u.y; M[6] = u.z; M[7] = 0;
	M[8] = -f.x; M[9] = -f.y; M[10] = -f.z; M[11] = 0;
	M[12] = 0; M[13] = 0; M[14] = 0; M[15] = 1;
	M[0] = s.x; M[4] = s.y; M[8] = s.z; M[12] = 0;
	M[1] = u.x; M[5] = u.y; M[9] = u.z; M[13] = 0;
	M[2] = -f.x; M[6] = -f.y; M[10] = -f.z; M[14] = 0;
	M[3] = 0; M[7] = 0; M[11] = 0; M[15] = 1;

	return	M;
}


void Point::CrossProduct(const Point p1, const Point p2, Point &p3)
{
	p3.x = p1.y * p2.z - p1.z * p2.y;
	p3.y = p1.z * p2.x - p1.x * p2.z;
	p3.z = p1.x * p2.y - p1.y * p2.x;
}


void Parser::Read(string filename)
{
	string line;
	m_file.open("level.nff");
	if(!m_file.is_open())
	{
		cout << "File not opened." << endl;
		exit(-1);
	}
	while(! m_file.eof())
	{
		string sCommand;
	//	getline(m_parser.m_file, line);
	//	cout << line << endl;
		m_file >> sCommand;
	//	cout << sCommand << endl;
		if(sCommand == "pp")
		{	cout << "It does not support pp" << endl; exit(-4);	}
		switch( sCommand[0] )
		{
			case ' ':            /* white space */
        		case '\t':
        		case '\n':
        		case '\f':
        		case '\r':
        		        continue;
			case '#':            /* comment */
				DoComment();
        		        break;
        		case 'v':            /* view point */
				DoViewPoint();
        		        break;
        		case 'l':            /* light source */
				DoLightSource();
        		        break;
        		case 'b':            /* background color */
				DoBackground();
        		        break;
        		case 'f':            /* fill material */
				DoFill();
        		        break;
        		case 'c':            /* cylinder or cone */
				DoCone();
        		        break;
        		case 's':            /* sphere */
				DoSphere();
        		        break;
			case 'p':            /* polygon or patch */
				DoPoly();
        		        break;
			case '\0' :
				break;
        		default:            /* unknown */
				cout << "unknown NFF primitive code" << endl;
				exit(1);
		}
	}
	g_view.SetView(m_view);
}

void Parser::DoViewPoint()
{
	string sCameraPos;
	m_file >> sCameraPos;
	if(sCameraPos != "from")
	{
		cout << "Token error >  Token != From. However, it is : " << sCameraPos << endl;
		exit(2);
	}
	//cout << sCameraPos << endl;
/*
        Point3   camera;
        Point3   lookAt;
        Vector3  cameraUp;
        float   angle;
        float   dNear;
        int   resolution;
*/
	m_file >> m_view.camera.x;	// From == Camera position
	//cout << m_view.camera.x << endl;
	m_file >> m_view.camera.y;
	//cout << m_view.camera.y << endl;
	m_file >> m_view.camera.z;
	//cout << m_view.camera.z << endl;

	string sLookAt;
	m_file >> sLookAt;
	if(sLookAt != "at")
	{
		cout << "The token must be \"at\", but it is " << sLookAt << endl;
		exit(3);
	}
	m_file >> m_view.lookAt.x;	// At == Camera looks at (look at position)
	m_file >> m_view.lookAt.y;
	m_file >> m_view.lookAt.z;

	string sCameraUp;
	m_file >> sCameraUp;
	if(sCameraUp != "up")
	{
		cout << "The token must be \"up\", but it is " << sCameraUp << endl;
		exit(4);
	}
	m_file >> m_view.cameraUp.x;	// Camera Up Vector Setting
	m_file >> m_view.cameraUp.y;
	m_file >> m_view.cameraUp.z;

	string sAngle;
	m_file >> sAngle;
	if(sAngle != "angle")
	{
		cout << "The token must be \"angle\", but it is " << sAngle << endl;
		exit(5);
	}
	m_file >> m_view.angle;		// Camera Angle.

	string sNear;
	m_file >> sNear;
	if(sNear != "hither")
	{
		cout << "The token must be \"hither\", but it is " << sNear << endl;
		exit(5);
	}
	m_file >> m_view.dNear;	// Distance to Near Plane Setting

	string sResolution;
	m_file >> sResolution;
	if(sResolution != "resolution")
	{
		cout << "The token must be \"resolution\", but it is " << sResolution << endl;
		exit(6);
	}
	m_file >> m_view.resolution.x;	// Resolution Setting
	m_file >> m_view.resolution.y;
}

void Parser::DoLightSource()
{
	m_file >> m_light.x;
	m_file >> m_light.y;
	m_file >> m_light.z; 
}

void Parser::DoBackground()
{
	m_file >> m_background.R;
	m_file >> m_background.G;
	m_file >> m_background.B;
}

void Parser::DoFill()
{
/*
	Color3	color;
	float	Kd;
	float	Ks;
	float	PhongPow;
	float	Transmission;
	float	indexRefraction;
*/
	m_file >> m_shading.color.R;
	m_file >> m_shading.color.G;
	m_file >> m_shading.color.B;
	m_file >> m_shading.Kd;
	m_file >> m_shading.Ks;
	m_file >> m_shading.PhongPow;
	m_file >> m_shading.Transmission;
	m_file >> m_shading.indexRefraction;
}

void Parser::DoCone()
{
}

void Parser::DoSphere()
{
	// # of spheres increments
/*	m_nSphere ++;
	// Now allocate memory for a sphere. But before that, parse it to extract data first.
	SPHERE sp;
	m_file >> sp.x;
	m_file >> sp.y;
	m_file >> sp.z;
	m_file >> sp.radius;

	m_lSphere.push_front(sp);*/
}

void Parser::DoPoly()
{
	int nVertices;
	m_file >> nVertices;
	if(nVertices < 3)
	{
		cout << "Polygons have more than 3 vertices." << endl;
		exit(-6);
	}
//	I modified the code because in this code I do not need triangles. I just get rectangles.
	if(nVertices == 4)
	{
		float vertices[4][3];
		TRI_STR tri;
		for(int i = 0 ; i < 4 ; i++)
		{
			for(int j = 0 ; j < 3 ; j++)
			{
				m_file >> vertices[i][j];
			}
		}
		for(int i = 0 ; i < 3 ; i++)
		{
			tri.a[i] = vertices[0][i];
		}
		for(int i = 0 ; i < 3 ; i++)
		{
			tri.b[i] = vertices[1][i];
		}
		for(int i = 0 ; i < 3 ; i++)
		{
			tri.c[i] = vertices[2][i];
		}
		for(int i = 0 ; i < 3 ; i++)
		{
			tri.d[i] = vertices[3][i];
		}
		m_lTri.push_front(tri);
	}

	 
/*	if(nVertices >= 3)
	{
		int nTriangles = nVertices - 3 + 1;
		float vertices[nVertices][3];
		for(int i = 0 ; i < nVertices ; i++)
		{
			m_file >> vertices[i][0];
			m_file >> vertices[i][1];
			m_file >> vertices[i][2];
		}
		int TriangleVertexNo[10][3] = { {0, 1, 2},
						{0, 2, 3},
						{0, 3, 4},
						{0, 4, 5},
						{0, 5, 6},
						{0, 6, 7},
						{0, 7, 8},
						{0, 8, 9},
						{0, 9, 10}};
		// Now I have all the vertex array.
		for(int i = 0 ; i < nTriangles ; i++)
		{
			// Get the 3 vertices for now.
			float a[3], b[3], c[3];
			a[0] = vertices[TriangleVertexNo[i][0]][0];
			a[1] = vertices[TriangleVertexNo[i][0]][1];
			a[2] = vertices[TriangleVertexNo[i][0]][2];
			b[0] = vertices[TriangleVertexNo[i][1]][0];
			b[1] = vertices[TriangleVertexNo[i][1]][1];
			b[2] = vertices[TriangleVertexNo[i][1]][2];
			c[0] = vertices[TriangleVertexNo[i][2]][0];
			c[1] = vertices[TriangleVertexNo[i][2]][1];
			c[2] = vertices[TriangleVertexNo[i][2]][2];
			// Now we have 3 vertices of the current triangle.

			// Time to add the triangle.
			m_nTetra ++;
			TETRA te;
			te.a[0] = a[0];
			te.a[1] = a[1];
			te.a[2] = a[2];
			te.b[0] = b[0];
			te.b[1] = b[1];
			te.b[2] = b[2];
			te.c[0] = c[0];
			te.c[1] = c[1];
			te.c[2] = c[2];

		// Compute N and also D as in Ax + By + Cz + D = 0
		// I compute D because I think I can use the equation to process polygons with 4/5 vertices.
			float b_a[3], c_a[3];
			MatrixSubtract1X3(te.b, te.a, b_a);
			MatrixSubtract1X3(te.c, te.a, c_a);
			CrossProduct(b_a, c_a, te.N);
			CrossProduct(te.N, te.a, te.D);
			te.D[0] = -te.D[0];
			te.D[1] = -te.D[1];
			te.D[2] = -te.D[2];

			m_lTetra.push_front(te);
		}
	}*/
}


void Parser::DoComment()
{
	char sComment[500];
	m_file.getline(sComment, 500);
}

void draw(void)
{
    /* clear old screen contents */
    glClearColor(0.5, 0.7, 0.9, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /* draw something */
	
	//turn on color and vertex drawing
	glEnableClientState(GL_COLOR_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	
//	glColorPointer(3, GL_FLOAT, 0, colors);
//	glVertexPointer(3,GL_FLOAT,0,verts);
	
	// draw the triangles
//    glDrawElements(GL_TRIANGLES,3,GL_UNSIGNED_BYTE,indices);
	
	//disable vertex and color arrays

	for(list<TRI_STR>::iterator it = g_parser.m_lTri.begin() ; it != g_parser.m_lTri.end() ; it++)
	{	// Color??? ==> m_shading.color.R, G, B
		TRI_STR tri = *it;
		glBegin(GL_QUADS);
			glColor3f(g_parser.m_shading.color.R, g_parser.m_shading.color.G, g_parser.m_shading.color.B);
			glVertex3f(tri.a[0], tri.a[1], tri.a[2]);
glColor3f(g_parser.m_shading.color.R - .3, g_parser.m_shading.color.G - 0.1, g_parser.m_shading.color.B + 0.3);
			glVertex3f(tri.b[0], tri.b[1], tri.b[2]);
glColor3f(g_parser.m_shading.color.R - .5, g_parser.m_shading.color.G - 0.3, g_parser.m_shading.color.B + 0.6);
			glVertex3f(tri.c[0], tri.c[1], tri.c[2]);
glColor3f(g_parser.m_shading.color.R - .6, g_parser.m_shading.color.G - 0.4, g_parser.m_shading.color.B + 0.9);
			glVertex3f(tri.d[0], tri.d[1], tri.d[2]);
		glEnd();
	}

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
    glutSwapBuffers();
}

void key(unsigned char k, int x, int y)
{
    switch (k)
    {
        case 27:			/* Escape: exit */
	       exit(0);

        case 'z':
            // z was pressed
            break;

        case 'x':
            // x was pressed
            break;
    }
}
