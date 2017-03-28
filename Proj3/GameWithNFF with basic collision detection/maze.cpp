/*
I got some of the code from the link and modified it.
http://nehe.gamedev.net/data/lessons/lesson.asp?lesson=10

I also referred to Wesley Griffin's code at <http://userpages.umbc.edu/~griffin5/cs435/01_cube.zip>

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
#define DIST_FROM_WALL	0.7

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
				cout << "hit ";
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

// loads the world from a text file.
void SetupWorld(string nff_filename) 
{
	string filename = nff_filename;
	g_parser.Read(nff_filename);
	// Now I need to display all the triangles.

 /*   float x, y, z, u, v;
    int vert;
    GLuint numtriangles;
    FILE *filein;        // file to load the world from
    char oneline[255];

    filein = fopen("Data/world.txt", "rt");

    readstr(filein, oneline);
    sscanf(oneline, "NUMPOLLIES %d\n", &numtriangles);

    sector1.numtriangles = numtriangles;
    sector1.triangle = (TRIANGLE *) malloc(sizeof(TRIANGLE)*numtriangles);
    
    for (loop = 0; loop < numtriangles; loop++) {
	for (vert = 0; vert < 3; vert++) {
	    readstr(filein,oneline);
	    sscanf(oneline, "%f %f %f %f %f", &x, &y, &z, &u, &v);
	    sector1.triangle[loop].vertex[vert].x = x;
	    sector1.triangle[loop].vertex[vert].y = y;
	    sector1.triangle[loop].vertex[vert].z = z;
	    sector1.triangle[loop].vertex[vert].u = u;
	    sector1.triangle[loop].vertex[vert].v = v;
	}
    }

    fclose(filein);*/
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

int nJumping = 0;
float jumpHeight = 0;
float jumpDelta = 0.1;
float jumpAcceleration = -0.0005;
/* The main drawing function. */
GLvoid DrawGLScene(GLvoid)
{
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
    //gluLookAt(0, 0, 0, 0, 0, 1, 0, 1, 0);
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
 /*   GLfloat xtrans, ztrans, ytrans;
    GLfloat sceneroty;

    // calculate translations and rotations.

    xtrans = -xpos;
    ytrans = 0.0;
    ztrans = -zpos;
    sceneroty = 0.0f - yrot;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		// Clear The Screen And The Depth Buffer
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    g_view.myOwnGluLookAt();
//gluLookAt(0, 0, 0, 0, 0, 1, 0, 1, 0);
//     gluLookAt(g_parser.m_view.camera.x, g_parser.m_view.camera.y, g_parser.m_view.camera.z,
//		g_parser.m_view.lookAt.x, g_parser.m_view.lookAt.y, g_parser.m_view.lookAt.z,
//		g_parser.m_view.cameraUp.x, g_parser.m_view.cameraUp.y, g_parser.m_view.cameraUp.z);

    if(nJumping == 1)
    {
	jumpDelta += jumpAcceleration;
	jumpHeight += jumpDelta/10;
	if(jumpHeight < 0)
	{
		nJumping = 0;
		jumpHeight = 0;
		jumpDelta = 0.1;
	}
	glutPostRedisplay();
    }
// Hacky +1 --> y
    g_view.myOwnGlRotate(g_AngleX, 1.0, 0.0, 0.0);
    g_view.myOwnGlRotate(g_AngleY, 0.0, 1.0, 0.0);
    g_view.myOwnGlRotate(g_AngleZ, 0.0, 0.0, 1.0);
    g_view.myOwnGlTranslate(xtrans, ytrans + jumpHeight + 1, ztrans);*/

//    glBindTexture(GL_TEXTURE_2D, texture[filter]);    // pick the texture.
//    numtriangles = sector1.numtriangles;

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
    float heading = 0.1;
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
	dir[0] = Vector(sin(yrot*piover180), 0.0, cos(yrot*piover180));	// forward
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
	hit_closest = collide.GetBarycentricCoordinate(ray_origin, dir[0]);
//	hit_closest = collide.GetClosestAmongDirections(ray_origin, dir, &nDir);
	if(hit_closest.closest_t < DIST_FROM_WALL)
		//Don't move!
		break;
//	hit_closest = collide.GetBarycentricCoordinate(ray_origin,  dir[0]); 
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
		//Don't move!
		break;
	xpos -= (float)sin(yrot*piover180) * heading;
	zpos -= (float)cos(yrot*piover180) * heading;
	break;

    case GLUT_KEY_LEFT: // look left
	//hit_closest = collide.GetClosestAmongDirections(ray_origin, dir, &nDir);
	hit_closest = collide.GetBarycentricCoordinate(ray_origin, dir[3]);
	if(hit_closest.closest_t < DIST_FROM_WALL)
		//Don't move!
		break;
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
		//Don't move!
		break;
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

int main(int argc, char **argv) 
{  
    /* load our world from disk */
// 	if(argc == 1)
// 	{
// 		cout << "You got to give me a nff file" << endl;
// 		exit(-6);
// 	}
    string nff_filename = "level.nff";
    SetupWorld(nff_filename);

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

			m_lTetra.push_front(te);
		}
	}
}


void Parser::DoComment()
{
	char sComment[500];
	m_file.getline(sComment, 500);
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

