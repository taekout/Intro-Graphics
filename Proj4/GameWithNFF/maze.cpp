/*
I got some of the code from the link and modified it.
http://nehe.gamedev.net/data/lessons/lesson.asp?lesson=10

I also referred to Wesley Griffin's code at <http://userpages.umbc.edu/~griffin5/cs435/01_cube.zip>

I was not sure about the shape of fractal, so I inquired Dr. Rheingans and she said it is okay.

*/
#include <GL/glut.h>    // Header File For The GLUT Parser 
#include <GL/gl.h>	// Header File For The OpenGL32 Parser
#include <GL/glu.h>	// Header File For The GLu32 Parser
#include <unistd.h>     // Header file for sleeping.
#include <stdio.h>      // Header file for standard file i/o.
#include <stdlib.h>     // Header file for malloc/free.
#include <math.h>       // Header file for trigonometric functions.
#include <iostream>
#include <list>
#include <fstream>
#include <time.h>
#include <sys/time.h>

#define SNOWFLAKES

using namespace std;

/* ascii codes for various special keys */
#define ESCAPE 27
#define PAGE_UP 73
#define PAGE_DOWN 81
#define UP_ARROW 72
#define DOWN_ARROW 80
#define LEFT_ARROW 75
#define RIGHT_ARROW 77

#define ZERO  1.0e-8
#define SPACE_END 100
#define HACKING_Y -1
#define DIST_FROM_WALL	1

// float g_m[16] = { 1.0, 0.0, 0.0, 0.0,
// 		0.0, 1.0, 0.0, 0.0,
// 		0.0, 0.0, 1.0, 0.0,
// 		0.0, 0.0, 0.0, 1.0};

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

typedef struct PLANE
{
	float a[3];
	float b[3];
	float c[3];
	float d[3];
	float N[3];
	float D[3];
} PLANE;

typedef struct PLANE5
{
	float a[3];
	float b[3];
	float c[3];
	float d[3];
	float e[3];
	float N[3];
	float D[3];
} PLANE5;

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


class Parser
{
public:
	ifstream	m_file;

	int m_nTetra;
	list<TETRA>	m_lTetra;
	Point		m_light;

	VIEW		m_view;
	COLOR3		m_background;
	SHADING		m_shading;

	Point		m_cameraVector;

	// daisy flower object circles
	virtual void Read(string filename);

	// Member Functions
	void DoViewPoint();
	void DoLightSource();
	void DoBackground();
	void DoFill();
	virtual void DoCone();
	virtual void DoSphere();
	virtual void DoPoly();
	void DoComment();
};

class ObjParser : public Parser
{
public:
/*
DAISY - P(3 vertices), Circle
CUBE - Polygon(4 vertices)
dodecahedron - P(5 v)
icosahedron - P(3 v)
octahedron - P(3 v)
tetrahedron - P(3 v)
swisscheese - P(3 v, 4 v), Circle(many)
*/
	// list of objects
	list<TETRA> m_lDaisyTetra;
	list<SPHERE> m_lDaisyCircle;
	list<PLANE> m_lCubePlane;
	list<PLANE5> m_lDodePlane5;
	list<TETRA> m_lIcosaTetra;
	list<TETRA> m_lOctaTetra;
	list<TETRA> m_lTetrahedronTetra;
	list<TETRA> m_lSwissTetra;
	list<PLANE> m_lSwissPlane;
	list<SPHERE> m_lSwissCircle;
	list<PLANE> m_lSnowflakePlane;
	list<Cylinder> m_lSnowflakeCylinder;
	list<PLANE> m_lSeashellPlane;
	list<TETRA> m_lPorcupineTetra;

	// Saved random Position of Objects
	int nDaisyObj;
	Point m_daisyPosition[30];
	int nCubeObj;
	Point m_cubePosition[30];
	int nDodeObj;
	Point m_DodePosition[30];
	int nIcosaObj;
	Point m_icosaPosition[30];
	int nOctaObj;
	Point m_octaPosition[30];
	int nTetrahedronObj;
	Point m_tetrahedronPosition[30];
	int nSwissObj;
	Point m_swissPosition[30];
	int nSnowflakeObj;
	Point m_snowflakePosition[30];
	int nSeashellObj;
	Point m_seashellPosition[30];
	int nPorcupineObj;
	Point m_porcupinePosition[30];

	ObjParser() {
		nDaisyObj = 0;
		nCubeObj = 0;
		nDodeObj = 0;
		nIcosaObj = 0;
		nOctaObj = 0;
		nTetrahedronObj = 0;
		nSwissObj = 0;
		nSnowflakeObj = 0;
		nSeashellObj = 0;
		nPorcupineObj = 0;
	}
	~ObjParser() {}

	// Later add snowflake
	// Parsing Function for each obj.
	void ReadObj(char filename[]);


	/*void ReadDaisyNFF();
	void ReadCubeNFF();
	void ReadDodecahedronNFF();
	void ReadIcosahedronNFF();
	void ReadOctahedronNFF();
	void ReadTetrahedronNFF();
	void ReadSwissNFF();*/

	virtual void DoSphere(string filename);
	virtual void DoPoly(string filename);
	virtual void DoCone(string filename);

	// Draw functions for each obj.
	void DrawDaisy();
	void drawCircle(double x, double y, double z, double radius, float colorR,float colorG,float colorB);
/*
nKind == 1 ==> Tetrahedron
nKind == 2 ==> Cube
nKind == 3 ==> octahedron
nKind == 4 ==> dodecahedron
nKind == 5 ==> icosahedron
*/
	void DrawPlatonic(int nKind);
	void DrawSwisscheese();
	void DrawSnowflake();
	void DrawSeashell();
	void DrawPorcupine();
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
	float *ComputeViewMatrix(Point **cameraUp = NULL, Point **cameraRight = NULL, Point **cameraViewDir = NULL);
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
	void myOwnGluLookAt();
	void myOwnGlTranslate(float x, float y, float z);
	void myOwnGlRotate(float angle, float x, float y, float z);

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

typedef Point Vector;
Viewing g_view;
Parser g_parser;
ObjParser g_objParser;
int nNff = -1;
char nffNames[20][30];

typedef struct RAY_INTERSECT{
	float closest_t;
	Point p;
	int index;
} RAY_INTERSECT;

class Collide
{
public:
	Collide() {}
	~Collide() {}

	// Member Functions
	Vector *GetNormal(int nIndex);
	RAY_INTERSECT GetBarycentricCoordinate(Point ray_origin_p, Vector ray_direction_v)
	{
		RAY_INTERSECT	ray_inter;
		ray_inter.closest_t = SPACE_END;
		ray_inter.p = Point(-1000, -1000, -1000);
		ray_inter.index = -1;
		float ray_origin[3] = {ray_origin_p.x, ray_origin_p.y, ray_origin_p.z};
		float ray_direction[3] = {ray_direction_v.x, ray_direction_v.y, ray_direction_v.z};
		float closest_t = SPACE_END;
		int index = 0;
		for(list<TETRA>::iterator it = g_parser.m_lTetra.begin() ; it != g_parser.m_lTetra.end() ; it++, index++)
		{
			// Compute A
			float a, b, c, d, e, f, g, h, i;
			a = (*it).a[0] - (*it).b[0];
			b = (*it).a[1] - (*it).b[1];
			c = (*it).a[2] - (*it).b[2];
			d = (*it).a[0] - (*it).c[0];
			e = (*it).a[1] - (*it).c[1];
			f = (*it).a[2] - (*it).c[2];
			g = ray_direction[0];
			h = ray_direction[1];
			i = ray_direction[2];
			float j, k, l;
			j = (*it).a[0] - ray_origin[0];
			k = (*it).a[1] - ray_origin[1];
			l = (*it).a[2] - ray_origin[2];

			// Compute M
			float M = a * (e * i - h *f) + b * (g * f - d * i) + c * (d * h - e * g);
			// if M == 0 ==> no solution ==> not intersect
			if(M == 0)
			{
//				cout << "The triangle is degenerate or parallel to the ray." << endl;
				continue;
				//return	-1;
			}
			// if M != 0 ==> get alpha, beta, t
			float t = f * (a * k - j * b) + e * (j * c - a * l) + d * (b * l - k * c);
			t /= M;
			t = -t;
			float beta = j * (e * i - h * f) + k * (g * f - d * i) + l * (d * h - e * g);
			beta /= M;
			float gamma = i * (a * k - j * b) + h * (j * c - a * l) + g * (b * l - k * c);
			gamma /= M;
			// Check they intersect inside the triangle.
			// check t  , gamma, beta ranges to see if it is really inside
			// If not inside, just continue; to let go of it and test some other triangles.
			if(gamma < 0 || gamma > 1)
				continue;
			if(beta < 0 || beta + gamma > 1)
				continue;
			if(t < 0)	// If t is behind the camera => does not intersect!
				continue;
			if(closest_t > t)
			{
				closest_t = t;
				ray_inter.closest_t = t;
				Point tmp(ray_direction_v.x * t, ray_direction_v.y * t, ray_direction_v.z * t); 
				ray_inter.p = ray_origin_p + tmp;
				ray_inter.index = index;
				//cout << "hit ";
			}
		}
		return	ray_inter;
	}
	// nDirection ==> 0 - forward, 1 - right, 2 - backward, 3 - left
	RAY_INTERSECT GetClosestAmongDirections(Point ray_origin, Vector v[4], int *nDirection);
	/*int TestIntersionPlane(const PLANE& plane,const Vector& position,const Vector& direction, double& lamda, Vector& pNormal)
	{
//		double DotProduct=direction.dot(plane._Normal);
		Point planeNormal(plane.N[0], plane.N[1], plane.N[2]); 
		double DotProduct = Point::DotProduct(direction, planeNormal);
		double l2;
	
		//determine if ray paralle to plane
		if( (DotProduct<ZERO)&&(DotProduct>-ZERO)   ) 
			return 0;
		l2=(plane._Normal.dot(plane._Position-position))/DotProduct;

		if (l2<-ZERO) 
			return 0;

		pNormal=plane._Normal;
		lamda=l2;
		return 1;
	}*/
};

/* The number of our GLUT window */
int window; 

GLuint loop;             // general loop variable
GLuint texture[3];       // storage for 3 textures;

int light = 0;           // lighting on/off
int blend = 0;        // blending on/off

GLfloat xrot;            // x rotation
GLfloat yrot;            // y rotation
GLfloat xspeed;          // x rotation speed
GLfloat yspeed;          // y rotation speed


//GLfloat lookupdown = 0.0; // Instead , I use GLfloat g_AngleX, g_AngleY, g_AngleZ;
GLfloat g_AngleX, g_AngleY, g_AngleZ;
int g_ButtonDown, g_LastX = -1, g_LastY = -1;
const float piover180 = 0.0174532925f;

float xpos, zpos, ypos;

//GLfloat camx = 0, camy = 0, camz = 0; // camera location.
GLfloat therotate;

GLfloat z=0.0f;                       // depth into the screen.

GLfloat LightAmbient[]  = {0.5f, 0.5f, 0.5f, 1.0f}; 
GLfloat LightDiffuse[]  = {1.0f, 1.0f, 1.0f, 1.0f}; 
GLfloat LightPosition[] = {-1.0f, -1.0f, -1.0f, -1.0f};

GLuint filter = 0;       // texture filtering method to use (nearest, linear, linear + mipmaps)

typedef struct {         // vertex coordinates - 3d and texture
    GLfloat x, y, z;     // 3d coords.
    GLfloat u, v;        // texture coords.
} VERTEX;

typedef struct {         // triangle
    VERTEX vertex[3];    // 3 vertices array
} TRIANGLE;

typedef struct {         // sector of a 3d environment
    int numtriangles;    // number of triangles in the sector
    TRIANGLE* triangle;  // pointer to array of triangles.
} SECTOR;

SECTOR sector1;

/* Image type - contains height, width, and data */
typedef struct {
    unsigned long sizeX;
    unsigned long sizeY;
    char *data;
} Image;

void CrossProduct(float M1[3], float M2[3], float resultM[3]);
void MatrixSubtract1X3(float M1[3], float M2[3], float resultM[3]);
void NormalizeVector(float M[3], float resultM[3]);

// degrees to radians...2 PI radians = 360 degrees
float rad(float angle)
{
    return angle * piover180;
}

// helper for SetupWorld.  reads a file into a string until a nonblank, non-comment line
// is found ("/" at the start indicating a comment); assumes lines < 255 characters long.
void readstr(FILE *f, char *string)
{
    do {
	fgets(string, 255, f); // read the line
    } while ((string[0] == '/') || (string[0] == '\n'));
    return;
}

Point RandomizePosition(void)
{
	int nPlaneNo = -1;
	nPlaneNo = rand() % 6; // [0,5] = nPlane

	float dx, dz;
	Point objPos;
	switch(nPlaneNo)
	{
	case 0:
		dx = (float)(rand() % 15);
		dx /= 10;	// dx = [0,1.4] <= 1.9
		dx += 0.25;
		objPos.x = -23.0 + dx;
		dz = (float)(rand() % 90);
		dz /= 10;	// dz = [0,8.9]; <= 9.9
		dz += 1.0/2;
		objPos.z = 10.0 + dz;
		objPos.y = 0.0f;
		break;
	case 1:
		dx = (float)(rand() % 175);
		dx /= 10;	// dx = [0,17.4] <= 19.9
		dx += 2.5/2;
		objPos.x = -21.0 + dx;
		dz = (float)(rand() % 15);
		dz /= 10;	// dz = [0,1.4]; <= 1.9
		dz += 0.5/2;
		objPos.z = 10.0 + dz;
		objPos.y = 0.0f;
		break;
	case 2:
		dx = (float)(rand() % 15);
		dx /= 10;	// dx = [0,1.4] <= 1.9
		dx += 0.5/2;
		objPos.x = -13.0 + dx;
		dz = (float)(rand() % 90);
		dz /= 10;	// dz = [0,8.9]; <= 9.9
		dz += 0.5;
		objPos.z = 0.0 + dz;
		objPos.y = 0.0f;
		break;
	case 3:
		dx = (float)(rand() % 15);
		dx /= 10;	// dx = [0,1.4] <= 1.9
		dx += 0.5/2;
		objPos.x = -1.0 + dx;
		dz = (float)(rand() % 90);
		dz /= 10;	// dz = [0,8.9]; <= 9.9
		dz += 1.0 / 2;
		objPos.z = 2.0 + dz;
		objPos.y = 0.0f;
		break;
	case 4:
		dx = (float)(rand() % 15);
		dx /= 10;	// dx = [0,1.4] <= 1.9
		dx += 0.5 / 2;
		objPos.x = 9.0 + dx;
		dz = (float)(rand() % 90);
		dz /= 10;	// dz = [0,8.9]; <= 9.9
		dz += 1.0 / 2;
		objPos.z = 2.0 + dz;
		objPos.y = 0.0f;
		break;
	case 5:
		dx = (float)(rand() % 200);
		dx /= 10;	// dx = [0,19.9] <= 21.9
		dx += 2.0 / 2;
		objPos.x = -11.0 + dx;
		dz = (float)(rand() % 15);
		dz /= 10;	// dz = [0,1.4]; <= 1.9
		dz += 0.5 / 2;
		objPos.z = 0.0 + dz;
		objPos.y = 0.0f;
		break;
	}
	return objPos;
}

// loads the world from a text file.
void SetupWorld(string nff_filename, bool bObj = false) 
{
	string filename = nff_filename;
	g_parser.Read(nff_filename);
	// Now I need to display all the triangles.Collide

	// Now Read Obj files and Draw it later.
	if(bObj == true)
	{
		// Process Obj.nff files
		// Let me just define draw files for each obj first.
		for(int i = 0 ; i < nNff ; i++)
		{
			g_objParser.ReadObj(nffNames[i]);
		}
	}

	return;
}

/*
 * getint and getshort are help functions to load the bitmap byte by byte on 
 * SPARC platform (actually, just makes the thing work on platforms of either
 * endianness, not just Intel's little endian)
 */
static unsigned int getint(FILE *fp)
{
  int c, c1, c2, c3;

  // get 4 bytes
  c = getc(fp);  
  c1 = getc(fp);  
  c2 = getc(fp);  
  c3 = getc(fp);
  
  return ((unsigned int) c) +   
    (((unsigned int) c1) << 8) + 
    (((unsigned int) c2) << 16) +
    (((unsigned int) c3) << 24);
}

static unsigned int getshort(FILE *fp)
{
  int c, c1;
  
  //get 2 bytes
  c = getc(fp);  
  c1 = getc(fp);

  return ((unsigned int) c) + (((unsigned int) c1) << 8);
}

// quick and dirty bitmap loader...for 24 bit bitmaps with 1 plane only.  
// See http://www.dcs.ed.ac.uk/~mxr/gfx/2d/BMP.txt for more info.
int ImageLoad(char *filename, Image *image) 
{
    FILE *file;
    unsigned long size;                 // size of the image in bytes.
    unsigned long i;                    // standard counter.
    unsigned short int planes;          // number of planes in image (must be 1) 
    unsigned short int bpp;             // number of bits per pixel (must be 24)
    char temp;                          // used to convert bgr to rgb color.

    // make sure the file is there.
    if ((file = fopen(filename, "rb"))==NULL) {
      printf("File Not Found : %s\n",filename);
      return 0;
    }

    // seek through the bmp header, up to the width/height:
    fseek(file, 18, SEEK_CUR);

    // No 100% errorchecking anymore!!!

    // read the width
    image->sizeX = getint (file);
    printf("Width of %s: %lu\n", filename, image->sizeX);
    
    // read the height 
    image->sizeY = getint (file);
    printf("Height of %s: %lu\n", filename, image->sizeY);
    
    // calculate the size (assuming 24 bits or 3 bytes per pixel).
    size = image->sizeX * image->sizeY * 3;

    // read the planes
    planes = getshort(file);
    if (planes != 1) {
	printf("Planes from %s is not 1: %u\n", filename, planes);
	return 0;
    }

    // read the bpp
    bpp = getshort(file);
    if (bpp != 24) {
      printf("Bpp from %s is not 24: %u\n", filename, bpp);
      return 0;
    }
	
    // seek past the rest of the bitmap header.
    fseek(file, 24, SEEK_CUR);

    // read the data. 
    image->data = (char *) malloc(size);
    if (image->data == NULL) {
	printf("Error allocating memory for color-corrected image data");
	return 0;	
    }

    if ((i = fread(image->data, size, 1, file)) != 1) {
	printf("Error reading image data from %s.\n", filename);
	return 0;
    }

    for (i=0;i<size;i+=3) { // reverse all of the colors. (bgr -> rgb)
	temp = image->data[i];
	image->data[i] = image->data[i+2];
	image->data[i+2] = temp;
    }

    // we're done.
    return 1;
}

// Load Bitmaps And Convert To Textures
GLvoid LoadGLTextures(GLvoid) 
{	
    // Load Texture
    Image *image1;
    
    // allocate space for texture
    image1 = (Image *) malloc(sizeof(Image));
    if (image1 == NULL) {
	printf("Error allocating space for image");
	exit(0);
    }

    if (!ImageLoad("Data/wall_bricks.bmp", image1)) {
	exit(1);
    }        

    // Create Textures	
    glGenTextures(3, &texture[0]);

    // nearest filtered texture
    glBindTexture(GL_TEXTURE_2D, texture[0]);   // 2d texture (x and y size)
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST); // scale cheaply when image bigger than texture
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST); // scale cheaply when image smalled than texture
    glTexImage2D(GL_TEXTURE_2D, 0, 3, image1->sizeX, image1->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, image1->data);

    // linear filtered texture
    glBindTexture(GL_TEXTURE_2D, texture[1]);   // 2d texture (x and y size)
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR); // scale linearly when image bigger than texture
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR); // scale linearly when image smalled than texture
    glTexImage2D(GL_TEXTURE_2D, 0, 3, image1->sizeX, image1->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, image1->data);

    // mipmapped texture
    glBindTexture(GL_TEXTURE_2D, texture[2]);   // 2d texture (x and y size)
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR); // scale linearly when image bigger than texture
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR_MIPMAP_NEAREST); // scale mipmap when image smalled than texture
    gluBuild2DMipmaps(GL_TEXTURE_2D, 3, image1->sizeX, image1->sizeY, GL_RGB, GL_UNSIGNED_BYTE, image1->data);
};

GLfloat fRand1[100] ;
GLfloat fRand2[100] ;
GLfloat fRand3[100] ;
/* A general OpenGL initialization function.  Sets all of the initial parameters. */
GLvoid InitGL(GLsizei Width, GLsizei Height)	// We call this right after our OpenGL window is created.
	{
	//    LoadGLTextures();                           // load the textures.
	//    glEnable(GL_TEXTURE_2D);                    // Enable texture mapping.
	
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);          // Set the blending function for translucency (note off at init time)
	glClearColor(g_parser.m_background.R, g_parser.m_background.G, g_parser.m_background.B, 1.0f);	// This Will Clear The Background Color To Black
	glClearDepth(1.0);				// Enables Clearing Of The Depth Buffer
	glDepthFunc(GL_LESS);                       // type of depth test to do.
	glEnable(GL_DEPTH_TEST);                    // enables depth testing.
	glEnable(GL_COLOR_MATERIAL);
	glShadeModel(GL_SMOOTH);			// Enables Smooth Color Shading
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();				// Reset The Projection Matrix
	
	gluPerspective(g_parser.m_view.angle,(GLfloat)1.0f, g_parser.m_view.dNear, 100.0f);	// Calculate The Aspect Ratio Of The Window
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	// set up lights.
	glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);
	glLightfv(GL_LIGHT1, GL_POSITION, LightPosition);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT1);
	for(int i = 0 ; i < 100 ; i++)
	{
		fRand1[i] = rand() % 100;
		fRand1[i] /= 1000.0;
		fRand2[i] = rand() % 100;
		fRand2[i] /= 1000.0;
		fRand3[i] = rand() % 100;
		fRand3[i] /= 1000.0;
	}
}


/* The function called when our window is resized (which shouldn't happen, because we're fullscreen) */
GLvoid ReSizeGLScene(GLsizei Width, GLsizei Height)
{
    if (Height==0)				// Prevent A Divide By Zero If The Window Is Too Small
	Height=1;

    glViewport(0, 0, Width, Height);		// Reset The Current Viewport And Perspective Transformation

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluPerspective(g_parser.m_view.angle,(GLfloat)1.0f, g_parser.m_view.dNear, 100.0f);	// Calculate The Aspect Ratio Of The Window
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

// jump variables
int nJumping = 0;
float jumpHeight = 0;
float jumpDelta = 0.1;
float jumpAcceleration = -0.0001;
// time variables
struct timeval tim1, tim2;
int bTime1Set = 0;
int bTime2Set = 0;
float mtime = 0.5f;
/* The main drawing function. */
GLvoid DrawGLScene(GLvoid)
{
    if(bTime1Set == 0)
    {
	gettimeofday(&tim1, NULL);
	bTime1Set = 1;
    }
    else
    {
	// First time set
	gettimeofday(&tim2, NULL);
	bTime2Set = 1;
	bTime1Set = 0;
	long seconds , useconds;
	seconds = tim2.tv_sec - tim1.tv_sec;
	useconds = tim2.tv_usec - tim1.tv_usec;
	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    }
    GLfloat x_m, y_m, z_m, u_m, v_m;
    GLfloat xtrans, ztrans, ytrans;
    GLfloat sceneroty;
    GLuint numtriangles;

    // calculate translations and rotations.
    xtrans = -xpos;
    ztrans = -zpos;
    ytrans = HACKING_Y;
    sceneroty = 360.0f - yrot;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		// Clear The Screen And The Depth Buffer
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //gluLookAt(10, 15, 5, 0, 0, 0, 0, 1, 0);
    g_view.myOwnGluLookAt();
//	g_view.ComputeViewMatrix();
//	glLoadMatrixf(g_view.viewMatrix);

    //g_view.ComputeViewMatrix();
    //glMultMatrixf(g_view.viewMatrix);

    //glRotatef(lookupdown, 1.0f, 0, 0);
    //glRotatef(sceneroty, 0, 1.0f, 0);
    g_view.myOwnGlRotate(g_AngleX, 1.0, 0.0, 0.0);
    g_view.myOwnGlRotate(g_AngleY, 0.0, 1.0, 0.0);
    g_view.myOwnGlRotate(g_AngleZ, 0.0, 0.0, 1.0);
//cout << "sceneroty " << sceneroty << endl;

    if(nJumping == 1)
    {
	jumpDelta += jumpAcceleration;
	jumpHeight += jumpDelta/100;
	if(jumpHeight < 0)
	{
		nJumping = 0;
		jumpHeight = 0;
		jumpDelta = 0.1;
	}
	glutPostRedisplay();
    }
//    glTranslatef(xtrans, ytrans - jumpHeight - 1, ztrans);
g_view.myOwnGlTranslate(xtrans, ytrans - jumpHeight, ztrans);

	for(list<TETRA>::iterator it = g_parser.m_lTetra.begin() ; it != g_parser.m_lTetra.end() ; it++)
	{	// Color??? ==> m_shading.color.R, G, B
		TETRA tet = *it;
		glBegin(GL_TRIANGLES);
			glColor3f(g_parser.m_shading.color.R, g_parser.m_shading.color.G, g_parser.m_shading.color.B);
			glVertex3f(tet.a[0], tet.a[1], tet.a[2]);
glColor3f(g_parser.m_shading.color.R - .5, g_parser.m_shading.color.G - 0.3, g_parser.m_shading.color.B + 0.3);
			glVertex3f(tet.b[0], tet.b[1], tet.b[2]);
glColor3f(g_parser.m_shading.color.R - .8, g_parser.m_shading.color.G - 0.6, g_parser.m_shading.color.B + 0.8);
			glVertex3f(tet.c[0], tet.c[1], tet.c[2]);
		glEnd();
	}

	// Start drawing object models. Start with Daisy.
	for(int i = 0 ; i < g_objParser.nDaisyObj ; i++)
	{
		glPushMatrix();
		glTranslatef(g_objParser.m_daisyPosition[i].x, g_objParser.m_daisyPosition[i].y, g_objParser.m_daisyPosition[i].z);
		g_objParser.DrawDaisy();
		glPopMatrix();
	}
/*
nKind == 1 ==> Tetrahedron
nKind == 2 ==> Cube
nKind == 3 ==> octahedron
nKind == 4 ==> dodecahedron
nKind == 5 ==> icosahedron
*/
	glTranslatef(0, 0.2, 0);
	for(int i = 0 ; i < g_objParser.nTetrahedronObj ; i++)
	{	// draw tetrahedrons
		glPushMatrix();
		glTranslatef(g_objParser.m_tetrahedronPosition[i].x, g_objParser.m_tetrahedronPosition[i].y, g_objParser.m_tetrahedronPosition[i].z);
		g_objParser.DrawPlatonic(1);
		glPopMatrix();
	}
	for(int i = 0 ; i < g_objParser.nCubeObj ; i++)
	{	// draw cubes
		glPushMatrix();
		glTranslatef(g_objParser.m_cubePosition[i].x, g_objParser.m_cubePosition[i].y, g_objParser.m_cubePosition[i].z);
		g_objParser.DrawPlatonic(2);
		glPopMatrix();
	}
	for(int i = 0 ; i < g_objParser.nCubeObj ; i++)
	{	// draw octahedron
		glPushMatrix();
		glTranslatef(g_objParser.m_octaPosition[i].x, g_objParser.m_octaPosition[i].y, g_objParser.m_octaPosition[i].z);
		g_objParser.DrawPlatonic(3);
		glPopMatrix();
	}
	for(int i = 0 ; i < g_objParser.nCubeObj ; i++)
	{	// draw dodecahedron
		glPushMatrix();
		glTranslatef(g_objParser.m_DodePosition[i].x, g_objParser.m_DodePosition[i].y, g_objParser.m_DodePosition[i].z);
		g_objParser.DrawPlatonic(4);
		glPopMatrix();
	}
	for(int i = 0 ; i < g_objParser.nCubeObj ; i++)
	{	// draw icosahedron
		glPushMatrix();
		glTranslatef(g_objParser.m_icosaPosition[i].x, g_objParser.m_icosaPosition[i].y, g_objParser.m_icosaPosition[i].z);
		g_objParser.DrawPlatonic(5);
		glPopMatrix();
	}
	for(int i = 0 ; i < g_objParser.nSwissObj ; i++)
	{	// draw swiss cheese
		glPushMatrix();
		glTranslatef(g_objParser.m_swissPosition[i].x, g_objParser.m_swissPosition[i].y, g_objParser.m_swissPosition[i].z);
		g_objParser.DrawSwisscheese();
		glPopMatrix();
	}

#ifdef SNOWFLAKES
	for(int i = 0 ; i < g_objParser.nSnowflakeObj ; i++)
	{	// draw snowflakes
		glPushMatrix();
		glTranslatef(g_objParser.m_snowflakePosition[i].x, g_objParser.m_snowflakePosition[i].y, g_objParser.m_snowflakePosition[i].z);
		g_objParser.DrawSnowflake();
		glPopMatrix();
	}
#endif
	for(int i = 0 ; i < g_objParser.nSeashellObj ; i++)
	{
		glPushMatrix();
		//glRotatef(180, 0, 0, 1);
		//glRotatef(-90, 0, 1, 0);
		glRotatef(90, 0, 0, 1);
		glTranslatef(g_objParser.m_seashellPosition[i].y, -g_objParser.m_seashellPosition[i].x, g_objParser.m_seashellPosition[i].z);
		g_objParser.DrawSeashell();
		glPopMatrix();
	}

	for(int i = 0 ; i < g_objParser.nPorcupineObj ; i++)
	{
		glPushMatrix();
		glColor3f(0, 0.5, 0.5);
		glTranslatef(g_objParser.m_porcupinePosition[i].x, g_objParser.m_porcupinePosition[i].y, g_objParser.m_porcupinePosition[i].z);
		g_objParser.DrawPorcupine();
		glPopMatrix();
	}

	// since this is double buffered, swap the buffers to display what just got drawn.
	glutSwapBuffers();
}


/* The function called whenever a normal key is pressed. */
void keyPressed(unsigned char key, int x, int y) 
{
    /* avoid thrashing this procedure */
    usleep(100);

    switch (key) {    
    case ESCAPE: // kill everything.
	/* exit the program...normal termination. */
	exit(1);                   	
	break; // redundant.

    case 'b': 
    case 'B': // switch the blending
	printf("B/b pressed; blending is: %d\n", blend);
	blend = blend ? 0 : 1;              // switch the current value of blend, between 0 and 1.
	if (blend) {
	    glEnable(GL_BLEND);
	    glDisable(GL_DEPTH_TEST);
	} else {
	    glDisable(GL_BLEND);
	    glEnable(GL_DEPTH_TEST);
	}
	printf("Blending is now: %d\n", blend);
	break;

    case 'f': 
    case 'F': // switch the filter
	printf("F/f pressed; filter is: %d\n", filter);
	filter++;                           // switch the current value of filter, between 0/1/2;
	if (filter > 2) {
	    filter = 0;
	}
	printf("Filter is now: %d\n", filter);
	break;

    case 'l': 
    case 'L': // switch the lighting
	printf("L/l pressed; lighting is: %d\n", light);
	light = light ? 0 : 1;              // switch the current value of light, between 0 and 1.
	if (light) {
	    glEnable(GL_LIGHTING);
	} else {
	    glDisable(GL_LIGHTING);
	}
	printf("Lighting is now: %d\n", light);
	break;

    case 32 :
	cout << "JUMP!" << endl;
	// Jump Status = 1;
	nJumping = 1;
	break;
    default:
      printf ("Key %d pressed. No action there yet.\n", key);
      break;
    }
}

/* The function called whenever a normal key is pressed. */
void specialKeyPressed(int key, int x, int y)
{
    /* avoid thrashing this procedure */
    usleep(100);
    float heading = 0.01 * mtime; // 0.01
	// Get World Cor Vector
	g_parser.m_view.cameraUp;
	Point worldUp(g_parser.m_view.cameraUp.x, g_parser.m_view.cameraUp.y, g_parser.m_view.cameraUp.z);
	worldUp.NormalizeVector();
	Point r;
	Point viewDir(-g_view.m_cameraOppositeViewDir ->x, -g_view.m_cameraOppositeViewDir ->y, -g_view.m_cameraOppositeViewDir ->z);
	viewDir.NormalizeVector(); 
	Point::CrossProduct(viewDir, worldUp, r);
	r.NormalizeVector();// r = World right, worldUp = World up, viewDir = World viewDir = traditionally -Z axis  

	// Get Ray
	Vector ray_origin(xpos, -HACKING_Y, zpos);
	Vector dir[4];
	dir[0] = Vector(sin(yrot*piover180), 0.0, cos(yrot*piover180));		// forward
	dir[1] = Vector(-cos(yrot*piover180), 0.0, sin(yrot*piover180));	// right
	dir[2] = Vector(-sin(yrot*piover180), 0.0, -cos(yrot*piover180));	// backward
	dir[3] = Vector(cos(yrot*piover180), 0.0, -sin(yrot*piover180));	// left
	for(int i = 0 ; i < 4 ; i++)
		dir[i].NormalizeVector();

	Collide collide;
	RAY_INTERSECT hit_closest;
	hit_closest.p = Point(-1000, -1000, -1000);
	hit_closest.closest_t = SPACE_END;
	hit_closest.index = -1;

	int nDir = -1;

    switch (key) {    
    case GLUT_KEY_UP: // walk forward (bob head)
///////////////////////////////////////////? !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  I got to check more than 2 walls. (more than 2 intersections)
	hit_closest = collide.GetBarycentricCoordinate(ray_origin, dir[0]);
//	hit_closest = collide.GetClosestAmongDirections(ray_origin, dir, &nDir);
	if(hit_closest.closest_t < DIST_FROM_WALL)
	{	//Move along the wall.
		Vector *normal = collide.GetNormal(hit_closest.index);
		Vector N(normal ->x, normal ->y, normal ->z);
		Vector wallDir;
		Vector upV(g_view.UP ->x, g_view.UP ->y, g_view.UP ->z); 
		Vector::CrossProduct(N, upV, wallDir);
		wallDir.NormalizeVector();
		if(Vector::DotProduct(wallDir, dir[0]) < 0)
			wallDir = Vector(-wallDir.x, -wallDir.y, -wallDir.z);
		xpos += wallDir.x * heading;
		zpos += wallDir.z * heading;
		break;
	}
// 	cout << "t:" << hit_closest.closest_t << " index:"<< hit_closest.index << " at:(" << hit_closest.p.x << "," << hit_closest.p.y << "," << hit_closest.p.z << ")"<< endl;
/*cout << ray_origin.x << " " << ray_origin.y << " " << ray_origin.z << endl;
cout << dir[0].x << " " << dir[0].y << " " << dir[0].z << endl;*/
	xpos += (float)sin(yrot*piover180) * heading; // piover180 --> a coefficient for conversion from degree to radian.
	zpos += (float)cos(yrot*piover180) * heading;
	break;

    case GLUT_KEY_DOWN: // walk back (bob head)
	//hit_closest = collide.GetClosestAmongDirections(ray_origin, dir, &nDir);
	hit_closest = collide.GetBarycentricCoordinate(ray_origin, dir[2]);
	if(hit_closest.closest_t < DIST_FROM_WALL)
	{	//Move along the wall.
		Vector *normal = collide.GetNormal(hit_closest.index);
		Vector N(normal ->x, normal ->y, normal ->z);
		Vector wallDir;
		Vector upV(g_view.UP ->x, g_view.UP ->y, g_view.UP ->z); 
		Vector::CrossProduct(N, upV, wallDir);
		wallDir.NormalizeVector();
		if(Vector::DotProduct(wallDir, dir[2]) < 0)
			wallDir = Vector(-wallDir.x, -wallDir.y, -wallDir.z);
		xpos += wallDir.x * heading;
		zpos += wallDir.z * heading;
		break;
	}
	xpos -= (float)sin(yrot*piover180) * heading;
	zpos -= (float)cos(yrot*piover180) * heading;
	break;

    case GLUT_KEY_LEFT: // look left
	//hit_closest = collide.GetClosestAmongDirections(ray_origin, dir, &nDir);
	hit_closest = collide.GetBarycentricCoordinate(ray_origin, dir[3]);
	if(hit_closest.closest_t < DIST_FROM_WALL)
	{	//Move along the wall.
		Vector *normal = collide.GetNormal(hit_closest.index);
		Vector N(normal ->x, normal ->y, normal ->z);
		Vector wallDir;
		Vector upV(g_view.UP ->x, g_view.UP ->y, g_view.UP ->z); 
		Vector::CrossProduct(N, upV, wallDir);
		wallDir.NormalizeVector();
		if(Vector::DotProduct(wallDir, dir[3]) < 0)
			wallDir = Vector(-wallDir.x, -wallDir.y, -wallDir.z);
		xpos += wallDir.x * heading;
		zpos += wallDir.z * heading;
		break;
	}
	xpos += (float)cos(yrot*piover180) * heading; // piover180 --> a coefficient for conversion from degree to radian.
	zpos -= (float)sin(yrot*piover180) * heading;
	break;

    case GLUT_KEY_RIGHT: // look right
/*	xpos += ((float)cos(yrot*piover180) * (g_view.m_cameraRight ->x) +
			(float)sin(yrot*piover180) *(g_view.m_cameraRight ->z)) * heading;
	ypos += (g_view.m_cameraRight ->y) * heading;
	zpos += ((float)(-sin(yrot*piover180)) * (g_view.m_cameraRight ->x) +
			(float)(cos(yrot*piover180)) * (g_view.m_cameraRight ->z)) * heading;*/
	//hit_closest = collide.GetClosestAmongDirections(ray_origin, dir, &nDir);
	hit_closest = collide.GetBarycentricCoordinate(ray_origin, dir[1]);
	if(hit_closest.closest_t < DIST_FROM_WALL)
	{	//Move along the wall.
		Vector *normal = collide.GetNormal(hit_closest.index);
		Vector N(normal ->x, normal ->y, normal ->z);
		Vector wallDir;
		Vector upV(g_view.UP ->x, g_view.UP ->y, g_view.UP ->z); 
		Vector::CrossProduct(N, upV, wallDir);
		wallDir.NormalizeVector();
		if(Vector::DotProduct(wallDir, dir[1]) < 0)
			wallDir = Vector(-wallDir.x, -wallDir.y, -wallDir.z);
		xpos += wallDir.x * heading;
		zpos += wallDir.z * heading;
		break;
	}
	xpos -= (float)cos(yrot*piover180) * heading; // piover180 --> a coefficient for conversion from degree to radian.
	zpos += (float)sin(yrot*piover180) * heading;
	break;

    default:
	printf ("Special key %d pressed. No action there yet.\n", key);
	break;
    }
	glutPostRedisplay();
}

// Handle any mouse button events.
void mouse(int button, int state, int x, int y) {
    // Record the mouse position when a button is pressed.
    g_LastX = x;
    g_LastY = y;
    g_ButtonDown = button;
//    g_ShiftDown = (glutGetModifiers() == GLUT_ACTIVE_SHIFT);
}

// Handle any mouse movement events.
// I gotta change yrot according to the mouse movement.
void motion(int x, int y) {
    float deltaY = x - g_LastX;
    float deltaX = y - g_LastY;

    deltaX /= 10;
    deltaY /= 10;

    switch(g_ButtonDown) {
    case GLUT_LEFT_BUTTON:                  // Rotate in the x-y plane
        g_AngleY -= deltaY;
	yrot += deltaY;				// Add mouse trackball
        g_AngleX -= deltaX;
        break;
/*    case GLUT_RIGHT_BUTTON:                 // Rotate in the x-z plane
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
        break;*/
    }

    g_LastX = x;
    g_LastY = y;
    glutPostRedisplay();
}

int main(int argc, char *argv[]) 
{
	srand(time(NULL));
	string nff_filename = "level.nff";

	cout << argv[2] << endl;
	string strConvert = argv[2];
	nNff = atoi(strConvert.c_str());
	cout << nNff << endl;
	if(nNff != -1)
	{
		for(int i = 0 ; i < nNff ; i++)
		{
			//nffNames[i] = argv[3 + i];
			strcpy(nffNames[i], argv[3 + i]);
		}
	}

	SetupWorld(nff_filename, true);
	
	/* Initialize GLUT state - glut will take any command line arguments that pertain to it or 
	X Windows - look at its documentation at http://reality.sgi.com/mjk/spec3/spec3.html */  
	glutInit(&argc, argv);  
	
	/* Select type of Display mode:   
	Double buffer 
	RGBA color
	Depth buffer 
	Alpha blending */  
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ALPHA);
	
	/* get a 640 x 480 window */
	glutInitWindowSize(g_parser.m_view.resolution.x, g_parser.m_view.resolution.y);
	
	/* the window starts at the upper left corner of the screen */
	glutInitWindowPosition(0, 0);  
	
	/* Open a window */  
	window = glutCreateWindow("Computer Graphics Project 3");  
	
	/* Register the function to do all our OpenGL drawing. */
	glutDisplayFunc(&DrawGLScene);
	
	/* Even if there are no events, redraw our gl scene. */
	glutIdleFunc(&DrawGLScene); 
	
	/* Register the function called when our window is resized. */
	glutReshapeFunc(&ReSizeGLScene);
	
	/* Register the function called when the keyboard is pressed. */
	glutKeyboardFunc(&keyPressed);
	
	/* Register the function called when special keys (arrows, page down, etc) are pressed. */
	glutSpecialFunc(&specialKeyPressed);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	
	/* Initialize our window. */
	InitGL(640, 480);
	/* Go fullscreen.  This is as soon as possible. */
	//glutFullScreen();
	
	
	/* Start Event Processing Engine */  
	glutMainLoop();  
	
	return 1;
}

void ObjParser::drawCircle(double x, double y, double z, double radius, float colorR,float colorG,float colorB)
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

	if(cameraUp != NULL && cameraRight != NULL && cameraOppositeViewDir != NULL)
	{
		*cameraUp = new Point(u.x, u.y, u.z);
		*cameraRight = new Point(s.x, s.y, s.z);
		*cameraOppositeViewDir = new Point(-f.x, -f.y, -f.z);
	}

	// Check the order !!! -> Vertical or Horizontal?
	float * M = new float[16];
//	M[0] = s.x; M[1] = s.y; M[2] = s.z; M[3] = 0;
//	M[4] = u.x; M[5] = u.y; M[6] = u.z; M[7] = 0;
//	M[8] = -f.x; M[9] = -f.y; M[10] = -f.z; M[11] = 0;
//	M[12] = 0; M[13] = 0; M[14] = 0; M[15] = 1;
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
	m_file.open(filename.c_str());
	if(!m_file.is_open())
	{
		cout << "File not opened." << endl;
		exit(-1);
	}
	while(! m_file.eof())
	{
		string sCommand;
		m_file >> sCommand;
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
/*
        Point3   camera;
        Point3   lookAt;
        Vector3  cameraUp;
        float   angle;
        float   dNear;
        int   resolution;
*/
	m_file >> m_view.camera.x;	// From == Camera position
	m_file >> m_view.camera.y;
	m_file >> m_view.camera.z;

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
	/*if(nVertices == 4)
	{
		float vertices[4][3];
		PLANE tri;
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
		// Compute N and also D as in Ax + By + Cz + D = 0
		// I compute D because I think I can use the equation to process polygons with 4/5 vertices.
		Point p1(tri.a[0], tri.a[1], tri.a[2]);
		Point p2(tri.b[0], tri.b[1], tri.b[2]);
		Point p3(tri.c[0], tri.c[1], tri.c[2]);
		Vector v1 = p2 - p1;
		Vector v2 = p3 - p1;

		v1.NormalizeVector();
		v2.NormalizeVector();
		Point N; Point D;
		Vector::CrossProduct(v1, v2, N);
		Vector::CrossProduct(N, p1, D);
		D.x = -D.x; D.y = -D.y; D.z = -D.z;

		m_lTri.push_front(tri);
	}*/

	if(nVertices >= 3)
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
			NormalizeVector(te.N, te.N);

			m_lTetra.push_front(te);
		}
	}
}

void ObjParser::DoCone(string filename)
{
	Cylinder cylin;
	m_file >> cylin.base_x;
	m_file >> cylin.base_y;
	m_file >> cylin.base_z;
	m_file >> cylin.base_radius;

	m_file >> cylin.apex_x;
	m_file >> cylin.apex_y;
	m_file >> cylin.apex_z;
	m_file >> cylin.apex_radius;

	if(filename == "snowflake.nff")
	{
		m_lSnowflakeCylinder.push_back(cylin);
	}
}

void ObjParser::DoSphere(string filename)
{
//Swiss and Daisy

	float x, y, z, radius;
	m_file >> x;
	m_file >> y;
	m_file >> z;
	m_file >> radius;
	SPHERE circle;
	circle.x = x;
	circle.y = y;
	circle.z = z;
	circle.radius = radius;
	if(filename == "daisy.nff")
	{
		m_lDaisyCircle.push_front(circle);
	}
	else if(filename == "swisscheese.nff")
	{
		m_lSwissCircle.push_front(circle);
	}
}

void ObjParser::DoPoly(string filename)
{
	int nVertices;
	m_file >> nVertices;
	if(nVertices < 3)
	{
		cout << "Polygons have more than 3 vertices." << endl;
		exit(-6);
	}
//	I modified the code because in this code I do not need triangles. I just get rectangles.
	if(nVertices == 5)
	{
		float vertices[5][3];
		PLANE5 tri;
		for(int i = 0 ; i < 5 ; i++)
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
		for(int i = 0 ; i < 3 ; i++)
		{
			tri.e[i] = vertices[4][i];
		}
		// Compute N and also D as in Ax + By + Cz + D = 0
		// I compute D because I think I can use the equation to process polygons with 4/5 vertices.
		Point p1(tri.a[0], tri.a[1], tri.a[2]);
		Point p2(tri.b[0], tri.b[1], tri.b[2]);
		Point p3(tri.c[0], tri.c[1], tri.c[2]);
		Vector v1 = p2 - p1;
		Vector v2 = p3 - p1;

		v1.NormalizeVector();
		v2.NormalizeVector();
		Point N; Point D;
		Vector::CrossProduct(v1, v2, N);
		Vector::CrossProduct(N, p1, D);
		D.x = -D.x; D.y = -D.y; D.z = -D.z;

		if(filename == "dodecahedron.nff")
			m_lDodePlane5.push_back(tri);
	}
	else if(nVertices == 4)
	{
		float vertices[4][3];
		PLANE tri;
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
		// Compute N and also D as in Ax + By + Cz + D = 0
		// I compute D because I think I can use the equation to process polygons with 4/5 vertices.
		Point p1(tri.a[0], tri.a[1], tri.a[2]);
		Point p2(tri.b[0], tri.b[1], tri.b[2]);
		Point p3(tri.c[0], tri.c[1], tri.c[2]);
		Vector v1 = p2 - p1;
		Vector v2 = p3 - p1;

		v1.NormalizeVector();
		v2.NormalizeVector();
		Point N; Point D;
		Vector::CrossProduct(v1, v2, N);
		Vector::CrossProduct(N, p1, D);
		D.x = -D.x; D.y = -D.y; D.z = -D.z;

		if(filename == "cube.nff")
			m_lCubePlane.push_back(tri);
		else if(filename == "swisscheese.nff")
			m_lSwissPlane.push_back(tri);
		else if(filename == "snowflake.nff")
			m_lSnowflakePlane.push_back(tri);
		else if(filename == "seashell.nff")
			m_lSeashellPlane.push_back(tri);
	}
	else if(nVertices == 3)
	{	// Read all vertices
		float vertices[3][3];
		for(int i = 0 ; i < nVertices ; i++)
		{
			m_file >> vertices[i][0];
			m_file >> vertices[i][1];
			m_file >> vertices[i][2];
		}
		// Now we have 3 vertices.
		// Time to add the triangle.
		m_nTetra ++;
		TETRA te;
		te.a[0] = vertices[0][0];
		te.a[1] = vertices[0][1];
		te.a[2] = vertices[0][2];
		te.b[0] = vertices[1][0];
		te.b[1] = vertices[1][1];
		te.b[2] = vertices[1][2];
		te.c[0] = vertices[2][0];
		te.c[1] = vertices[2][1];
		te.c[2] = vertices[2][2];

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
		NormalizeVector(te.N, te.N);

		if(filename == "daisy.nff")
			m_lDaisyTetra.push_back(te);
		else if(filename == "icosahedron.nff")
			m_lIcosaTetra.push_back(te);
		else if(filename == "octahedron.nff")
			m_lOctaTetra.push_back(te);
		else if(filename == "tetrahedron.nff")
			m_lTetrahedronTetra.push_back(te);
		else if(filename == "swisscheese.nff")
			m_lSwissTetra.push_back(te);
		else if(filename == "porcupine.nff")
			m_lPorcupineTetra.push_back(te);
	}
}


void Parser::DoComment()
{
	char sComment[500];
	m_file.getline(sComment, 500);
}

void ObjParser::ReadObj(char filename[])
{
	m_file.open(filename);
	if(!m_file.is_open())
	{
		cout << "File not opened." << endl;
		exit(-1);
	}

	string objNffName;
	objNffName = filename;
	if(objNffName == "daisy.nff")
	{
		m_daisyPosition[nDaisyObj] = RandomizePosition();
		nDaisyObj ++;
	}
	else if(objNffName == "tetrahedron.nff")
	{
		m_tetrahedronPosition[nTetrahedronObj] = RandomizePosition();
		nTetrahedronObj ++;
	}
	else if(objNffName == "cube.nff")
	{
		m_cubePosition[nCubeObj] = RandomizePosition();
		nCubeObj ++;
	}
	else if(objNffName == "dodecahedron.nff")
	{
		m_DodePosition[nDodeObj] = RandomizePosition();
		nDodeObj ++;
	}
	else if(objNffName == "icosahedron.nff")
	{
		m_icosaPosition[nIcosaObj] = RandomizePosition();
		nIcosaObj ++;
	}
	else if(objNffName == "octahedron.nff")
	{
		m_octaPosition[nOctaObj] = RandomizePosition();
		nOctaObj ++;
	}
	else if(objNffName == "swisscheese.nff")
	{
		m_swissPosition[nSwissObj] = RandomizePosition();
		nSwissObj ++;
	}
	else if(objNffName == "snowflake.nff")
	{
		m_snowflakePosition[nSnowflakeObj] = RandomizePosition();
		nSnowflakeObj ++;
	}
	else if(objNffName == "seashell.nff")
	{
		m_seashellPosition[nSeashellObj] = RandomizePosition();
		nSeashellObj ++;
	}
	else if(objNffName == "porcupine.nff")
	{
		m_porcupinePosition[nPorcupineObj] = RandomizePosition();
		nPorcupineObj ++;
	}



	while(! m_file.eof())
	{
		string sCommand;
		m_file >> sCommand;
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
				DoCone(filename);
        		        break;
        		case 's':            /* sphere */
				DoSphere(filename);
        		        break;
			case 'p':            /* polygon or patch */
				DoPoly(filename);
        		        break;
			case '\0' :
				break;
        		default:            /* unknown */
				cout << "unknown NFF primitive code" << endl;
				cout << sCommand << endl;
				exit(1);
		}
	}
	m_file.close();
}

// Draw functions for each obj.
void ObjParser::DrawDaisy()
{
	glTranslatef(0, 0.4, 0);

	glColor3f(0, 1.0, 0);
	int index = 0;
	for(list<TETRA>::iterator it = m_lDaisyTetra.begin() ; it != m_lDaisyTetra.end() ; it++, index++)
	{
		if(index == 8)
			glColor3f(1.0, 0.0, 0.0);
		TETRA te = *it;
		glBegin(GL_TRIANGLES);
			glVertex3fv(te.a);
			glVertex3fv(te.b);
			glVertex3fv(te.c);
		glEnd();
	}
	for(list<SPHERE>::iterator it = m_lDaisyCircle.begin() ; it != m_lDaisyCircle.end() ; it++)
	{
		SPHERE sp;
		sp = *it;
		drawCircle(sp.x, sp.y, sp.z, sp.radius, 0.5, 0, 0.5);
	}
}

/*
nKind == 1 ==> Tetrahedron
nKind == 2 ==> Cube
nKind == 3 ==> octahedron
nKind == 4 ==> dodecahedron
nKind == 5 ==> icosahedron
*/
void ObjParser::DrawPlatonic(int nKind)
{

// 	list<TETRA> m_lDaisyTetra;
// 	list<SPHERE> m_lDaisyCircle;
// 	list<PLANE> m_lCubePlane;
// 	list<PLANE5> m_lDodePlane5;
// 	list<TETRA> m_lIcosaTetra;
// 	list<TETRA> m_lOctaTetra;
// 	list<TETRA> m_lTetrahedronTetra;
// 	list<TETRA> m_lSwissTetra;
// 	list<PLANE> m_lSwissPlane;
// 	list<SPHERE> m_lSwissCircle;


	float variance = 0.3;
	glColor3f(0.9, 0.5, 0);
	if(nKind == 1) // Tetrahedron
	{// 4 poly , 3 vertices per poly
		for(list<TETRA>::iterator it = m_lTetrahedronTetra.begin() ; it != m_lTetrahedronTetra.end() ; it++)
		{
			TETRA tet = *it;
			glBegin(GL_TRIANGLES);
				glVertex3fv(tet.a);
				glColor3f(1.0, 0.7, 0);
				glVertex3fv(tet.b);
				glColor3f(1.0 - variance * 1, 0.7 - variance * 0.3 * 1, 0 + variance * 1 *0.2);
				glVertex3fv(tet.c);
				glColor3f(1.0 - variance * 2, 0.7 - variance * 0.3 * 2, 0 + variance * 2 *0.2);
			glEnd();
		}
	}
	else if(nKind == 2)	// Cube
	{	// 6 poly, 4 vertices per poly
		glBegin(GL_QUADS);
		for(list<PLANE>::iterator it = m_lCubePlane.begin() ; it != m_lCubePlane.end() ; it++)
		{
			PLANE pl = *it;
			glVertex3fv(pl.a);
			glColor3f(1.0, 0.7, 0);
			glVertex3fv(pl.b);
			glColor3f(1.0 - variance * 1, 0.7 - variance * 0.3 * 1, 0 + variance * 1 *0.2);
			glVertex3fv(pl.c);
			glColor3f(1.0 - variance * 2, 0.7 - variance * 0.3 * 2, 0 + variance * 2 *0.2);
			glVertex3fv(pl.d);
			glColor3f(1.0 - variance * 3, 0.7 - variance * 0.3 * 3, 0 + variance * 3 *0.2);
		}
		glEnd();
	}
	else if(nKind == 3) // octahedron
	{	// 8 poly, 3 vertices per poly.
		glBegin(GL_TRIANGLES);
		for(list<TETRA>::iterator it = m_lOctaTetra.begin() ; it != m_lOctaTetra.end() ; it++)
		{
			TETRA tet = *it;
			glVertex3fv(tet.a);
			glColor3f(1.0, 0.7, 0);
			glVertex3fv(tet.b);
			glColor3f(1.0 - variance * 1, 0.7 - variance * 0.3 * 1, 0 + variance * 1 *0.2);
			glVertex3fv(tet.c);
		}
		glEnd();
	}
	else if(nKind == 4) // dodecahedron
	{	// 12 poly , 5 vertices per poly.
		for(list<PLANE5>::iterator it = m_lDodePlane5.begin() ; it != m_lDodePlane5.end() ; it++)
		{
			glBegin(GL_POLYGON);
			PLANE5 pl5 = *it;
			glVertex3fv(pl5.a);
			glColor3f(1.0, 0.7, 0);
			glVertex3fv(pl5.b);
			glColor3f(1.0 - variance * 1, 0.7 - variance * 0.3 * 1, 0 + variance * 1 *0.2);
			glVertex3fv(pl5.c);
			glColor3f(1.0 - variance * 2, 0.7 - variance * 0.3 * 2, 0 + variance * 2 *0.2);
			glVertex3fv(pl5.d);
			glColor3f(1.0 - variance * 3, 0.7 - variance * 0.3 * 3, 0 + variance * 3 *0.2);
			glVertex3fv(pl5.e);
			glColor3f(1.0 - variance * 4, 0.7 - variance * 0.3 * 4, 0 + variance * 4 *0.2);
			glEnd();
		}
	}
	else if(nKind == 5) //icosahedron
	{
		glBegin(GL_TRIANGLES);
		for(list<TETRA>::iterator it = m_lIcosaTetra.begin() ; it != m_lIcosaTetra.end() ; it++)
		{
			TETRA tet = *it;
			glVertex3fv(tet.a);
			glColor3f(1.0, 0.7, 0);
			glVertex3fv(tet.b);
			glColor3f(1.0 - variance * 1, 0.7 - variance * 0.3 * 1, 0 + variance * 1 *0.2);
			glVertex3fv(tet.c);
			glColor3f(1.0 - variance * 2, 0.7 - variance * 0.3 * 2, 0 + variance * 2 *0.2);
		}
		glEnd();
	}
}

/*
swisscheese color -> 0.17254902 0.654901961  0.917647059
hole - > 0.164705882 1 0.917647059
MyVector4 cheese[2][3];
MyVector4 cheeseSide[3][4];
MyVector4 hole[2][6];
float holeRadius = 0.2;
*/
void ObjParser::DrawSwisscheese()
{
//234 193 82
	glColor3f(234.0/255, 193.0/255, 82.0/255);

	for(list<PLANE>::iterator it = m_lSwissPlane.begin() ; it != m_lSwissPlane.end() ; it++)
	{
		PLANE pl = *it;
		glBegin(GL_QUADS);
			glVertex3fv(pl.a);
			glVertex3fv(pl.b);
			glVertex3fv(pl.c);
			glVertex3fv(pl.d);
		glEnd();
	}
	for(list<TETRA>::iterator it = m_lSwissTetra.begin() ; it != m_lSwissTetra.end() ; it++)
	{
		TETRA te = *it;
		glBegin(GL_TRIANGLES);
			glVertex3fv(te.a);
			glVertex3fv(te.b);
			glVertex3fv(te.c);
		glEnd();
	}
	for(list<SPHERE>::iterator it = m_lSwissCircle.begin() ; it != m_lSwissCircle.end() ; it++)
	{
		SPHERE sp;
		sp = *it;
		drawCircle(sp.x, sp.y, sp.z, sp.radius, 230.0/255,162.0/255, 0);
	}
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

void ObjParser::DrawPorcupine(void)
{
	glBegin(GL_TRIANGLES);
	for(list<TETRA>::iterator it = m_lPorcupineTetra.begin() ; it != m_lPorcupineTetra.end() ; it++)
	{
		TETRA te = *it;
		float *N;
		te.a[0] /= 10; te.a[1] /= 10; te.a[2] /= 10;
		te.b[0] /= 10; te.b[1] /= 10; te.b[2] /= 10;
		te.c[0] /= 10; te.c[1] /= 10; te.c[2] /= 10;
		Triangle tri_class(te.a , te.b, te.c);
		N = ComputeNormal(tri_class);
		glNormal3fv(N);
		glVertex3fv(te.a);
		glVertex3fv(te.b);
		glVertex3fv(te.c);
	}
	glEnd();
}

GLdouble red = 0.670588235;
GLdouble green = 0.784313725;
GLdouble blue = 0.949019608;
GLdouble dRed = (0.078431373 - 0.670588235) / 2200.0;
GLdouble dGreen = (0.121568627 - 0.784313725) / 2200.0;
GLdouble dBlue = (0.37254902 - 0.949019608) / 2200.0;
void ObjParser::DrawSeashell()
{
	red = 0.670588235;
	green = 0.784313725;
	blue = 0.949019608;
	int i = 0;
	for(list<PLANE>::iterator it = m_lSeashellPlane.begin() ; it != m_lSeashellPlane.end() ; it++, i++)
	{
//from :171 200 242  0.670588235 0.784313725 0.949019608
//to : 20 31 95 0.078431373 0.121568627 0.37254902
		if(i == 100) i = 0;
		PLANE pl = *it;
		glBegin(GL_QUADS);
			glColor3f(red + fRand1[i], green + fRand2[i], blue + fRand3[i]);
			glVertex3fv(pl.a);
			glVertex3fv(pl.b);
			glVertex3fv(pl.c);
			glVertex3fv(pl.d);
		glEnd();
		red += dRed;
		green += dGreen;
		blue += dBlue;
	}
}

void ObjParser::DrawSnowflake()
{
/*	for(list<Cylinder>::iterator it = m_lSnowflakeCylinder.begin() ; it != m_lSnowflakeCylinder.end() ; it++)
	{
		Cylinder cy = *it;
		glPushMatrix();
		glTranslatef(cy.base_x, cy.base_y, cy.base_z);//1
		glColor3f(0.7, 0.3, 0.4);
		GLUquadric *qobj = gluNewQuadric();
		gluCylinder(qobj, cy.base_radius, cy.apex_radius, 0.3, 24, 24);
		glPopMatrix();
	}*/

	for(list<PLANE>::iterator it = m_lSnowflakePlane.begin() ; it != m_lSnowflakePlane.end() ; it++)
	{
		PLANE pl = *it;
		glBegin(GL_QUADS);
			glColor3f(0.7, 0.3, 0.4);
			glVertex3fv(pl.a);
			glVertex3fv(pl.b);
//			glColor3f(0.7, 0.6, 0.5);
			glVertex3fv(pl.c);
			glVertex3fv(pl.d);
		glEnd();
	}
}

void CrossProduct(float M1[3], float M2[3], float resultM[3])
{
	resultM[0] = M1[1] * M2[2] - M1[2] * M2[1];
	resultM[1] = M1[2] * M2[0] - M1[0] * M2[2];
	resultM[2] = M1[0] * M2[1] - M1[1] * M2[0];
}

void MatrixSubtract1X3(float M1[3], float M2[3], float resultM[3])
{
	resultM[0] = M1[0] - M2[0];
	resultM[1] = M1[1] - M2[1];
	resultM[2] = M1[2] - M2[2];
}

void NormalizeVector(float M[3], float resultM[3])
{
	float length = M[0] * M[0] + M[1] * M[1] + M[2] * M[2];
	length = sqrt(length);
	resultM[0] = M[0] / length;
	resultM[1] = M[1] / length;
	resultM[2] = M[2] / length;
}

void Viewing::myOwnGluLookAt()
{
	g_view.ComputeViewMatrix();
	glLoadMatrixf(g_view.viewMatrix);
//	glTranslatef(-g_view.cameraPos ->x, -g_view.cameraPos ->y, -g_view.cameraPos ->z);
	myOwnGlTranslate(g_view.cameraPos ->x, g_view.cameraPos ->y, g_view.cameraPos ->z);
}

void Viewing::myOwnGlTranslate(float x, float y, float z)
{
	GLfloat *m = new GLfloat[16];
	m[0] = 1; m[4] = 0; m[8] = 0; m[12] = x;
	m[1] = 0; m[5] = 1; m[9] = 0; m[13] = y;
	m[2] = 0; m[6] = 0; m[10] = 1; m[14] = z;
	m[3] = 0; m[7] = 0; m[11] = 0; m[15] = 1;
	glMultMatrixf(m);
}

// I dont want to do 4X4 matrix multiplication. Man... 
void Viewing::myOwnGlRotate(float angle, float x, float y, float z)
{
	float m[16];
	float c = cos(angle * piover180);
	float s = sin(angle * piover180);
	Point axis(x, y, z);
	axis.NormalizeVector();
	x = axis.x; y = axis.y; z = axis.z;
	m[0] = x * x * (1-c) + c; m[4] = x * y * (1-c) - z *s ; m[8] = x*z *(1-c)+y * s  ; m[12] = 0;
	m[1] = y * x *(1-c)+ z*s; m[5] = y * y *(1-c) + c     ; m[9] = y*z * (1-c) - x*s ; m[13] = 0;
	m[2] = x*z*(1-c) - y*s  ; m[6] = y * z * (1-c)+ x*s   ; m[10] = z*z*(1-c) + c    ; m[14] = 0;
	m[3] = 0                ; m[7] = 0                    ; m[11] = 0                ; m[15] = 1;

	glMultMatrixf(m);

}

// I am not using this function for this program
RAY_INTERSECT Collide::GetClosestAmongDirections(Point ray_origin, Vector v[4], int *nDirection)
{
	RAY_INTERSECT rayinfo[4];
	RAY_INTERSECT returnInfo;
	returnInfo.closest_t = SPACE_END;
	returnInfo.index = -1;
	returnInfo.p = Point(-1000, -1000, -1000);
	for(int i = 0 ; i < 4 ; i++)
	{
		rayinfo[i] = GetBarycentricCoordinate(ray_origin, v[i]);
cout << "t:" << rayinfo[i].closest_t << " index:"<< rayinfo[i].index << " at:(" << rayinfo[i].p.x << "," << rayinfo[i].p.y << "," << rayinfo[i].p.z << ")"<< endl;
	}
	float closest_t = SPACE_END;
	for(int i = 0 ; i < 4 ; i++)
	{
		if(rayinfo[i].closest_t < closest_t)
		{
			closest_t = rayinfo[i].closest_t;
			*nDirection = i;
			returnInfo.p = rayinfo[i].p;
			returnInfo.closest_t = rayinfo[i].closest_t;
			returnInfo.index = rayinfo[i].index;
		}
	}
	return	returnInfo;
}

Vector *Collide::GetNormal(int nIndex)
{
	int index = 0;
	for(list<TETRA>::iterator it = g_parser.m_lTetra.begin() ; it != g_parser.m_lTetra.end() ; it++, index++)
	{
		if(index == nIndex)
		{
			TETRA tet= *it;
			Vector *v = new Vector(tet.N[0], tet.N[1], tet.N[2]);
			return v;
		}
	}
	return	NULL;
}




