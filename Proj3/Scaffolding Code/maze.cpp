/*
I got some of the code from the link and modified it.
http://nehe.gamedev.net/data/lessons/lesson.asp?lesson=10

I also referred to Wesley Griffin's code at <http://userpages.umbc.edu/~griffin5/cs435/01_cube.zip>

*/
#include <GL/glut.h>    // Header File For The GLUT Library 
#include <GL/gl.h>	// Header File For The OpenGL32 Library
#include <GL/glu.h>	// Header File For The GLu32 Library
#include <unistd.h>     // Header file for sleeping.
#include <stdio.h>      // Header file for standard file i/o.
#include <stdlib.h>     // Header file for malloc/free.
#include <math.h>       // Header file for trigonometric functions.
#include <iostream>

using namespace std; 

/* ascii codes for various special keys */
#define ESCAPE 27
#define PAGE_UP 73
#define PAGE_DOWN 81
#define UP_ARROW 72
#define DOWN_ARROW 80
#define LEFT_ARROW 75
#define RIGHT_ARROW 77

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

class Viewing{
public:
	Viewing()
	{
		centerPos = new Point(0, 0, 0);
		cameraPos = new Point(0, 0, -1);
		UP = new Point(0, 1, 0);
		viewMatrix = new float[16];	
	}
	~Viewing()
	{
		delete centerPos;
		delete cameraPos;
		delete viewMatrix;
		delete UP;
	}

// Member Variables
	Point *centerPos;
	Point *cameraPos;
	Point *UP;
	float *viewMatrix;
// Member Functions
	float *ComputeViewMatrix();
};



Viewing g_view;

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

GLfloat walkbias = 0;
GLfloat walkbiasangle = 0;

//GLfloat lookupdown = 0.0; // Instead , I use GLfloat g_AngleX, g_AngleY, g_AngleZ;
GLfloat g_AngleX, g_AngleY, g_AngleZ;
int g_ButtonDown, g_LastX = -1, g_LastY = -1;
const float piover180 = 0.0174532925f;

float heading, xpos, zpos;

//GLfloat camx = 0, camy = 0, camz = 0; // camera location.
GLfloat therotate;

GLfloat z=0.0f;                       // depth into the screen.

GLfloat LightAmbient[]  = {0.5f, 0.5f, 0.5f, 1.0f}; 
GLfloat LightDiffuse[]  = {1.0f, 1.0f, 1.0f, 1.0f}; 
GLfloat LightPosition[] = {0.0f, 0.0f, 2.0f, 1.0f};

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
void SetupWorld() 
{
    float x, y, z, u, v;
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

    fclose(filein);
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
    LoadGLTextures();                           // load the textures.
    glEnable(GL_TEXTURE_2D);                    // Enable texture mapping.

    glBlendFunc(GL_SRC_ALPHA, GL_ONE);          // Set the blending function for translucency (note off at init time)
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);	// This Will Clear The Background Color To Black
    glClearDepth(1.0);				// Enables Clearing Of The Depth Buffer
    glDepthFunc(GL_LESS);                       // type of depth test to do.
    glEnable(GL_DEPTH_TEST);                    // enables depth testing.
    glShadeModel(GL_SMOOTH);			// Enables Smooth Color Shading

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();				// Reset The Projection Matrix

    gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);	// Calculate The Aspect Ratio Of The Window

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // set up lights.
    glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);
    glLightfv(GL_LIGHT1, GL_POSITION, LightPosition);
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

    gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

int nJumping = 0;
float jumpHeight = 0;
float jumpDelta = 0.1;
float jumpAcceleration = -0.001;
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
    ytrans = -walkbias-0.25f;
    sceneroty = 360.0f - yrot;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		// Clear The Screen And The Depth Buffer
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0, 0, 0, 0, 0, 1, 0, 1, 0);
//	g_view.ComputeViewMatrix();
//	glLoadMatrixf(g_view.viewMatrix);

    //g_view.ComputeViewMatrix();
    //glMultMatrixf(g_view.viewMatrix);

    //glRotatef(lookupdown, 1.0f, 0, 0);
    //glRotatef(sceneroty, 0, 1.0f, 0);
    glRotatef(g_AngleX, 1.0, 0.0, 0.0);
    glRotatef(g_AngleY, 0.0, 1.0, 0.0);
    glRotatef(g_AngleZ, 0.0, 0.0, 1.0);
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
    glTranslatef(xtrans, ytrans - jumpHeight, ztrans);    

    glBindTexture(GL_TEXTURE_2D, texture[filter]);    // pick the texture.

    numtriangles = sector1.numtriangles;

    for (loop=0; loop<numtriangles; loop++) {        // loop through all the triangles
	glBegin(GL_TRIANGLES);		
	glNormal3f( 0.0f, 0.0f, 1.0f);
	
	x_m = sector1.triangle[loop].vertex[0].x;
	y_m = sector1.triangle[loop].vertex[0].y;
	z_m = sector1.triangle[loop].vertex[0].z;
	u_m = sector1.triangle[loop].vertex[0].u;
	v_m = sector1.triangle[loop].vertex[0].v;
	glTexCoord2f(u_m,v_m); 
	glVertex3f(x_m,y_m,z_m);
	
	x_m = sector1.triangle[loop].vertex[1].x;
	y_m = sector1.triangle[loop].vertex[1].y;
	z_m = sector1.triangle[loop].vertex[1].z;
	u_m = sector1.triangle[loop].vertex[1].u;
	v_m = sector1.triangle[loop].vertex[1].v;
	glTexCoord2f(u_m,v_m); 
	glVertex3f(x_m,y_m,z_m);
	
	x_m = sector1.triangle[loop].vertex[2].x;
	y_m = sector1.triangle[loop].vertex[2].y;
	z_m = sector1.triangle[loop].vertex[2].z;
	u_m = sector1.triangle[loop].vertex[2].u;
	v_m = sector1.triangle[loop].vertex[2].v;
	glTexCoord2f(u_m,v_m); 
	glVertex3f(x_m,y_m,z_m);	
	
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

    switch (key) {    
/*    case GLUT_KEY_PAGE_UP: // tilt up
	z -= 0.2f;
	lookupdown -= 1.0f;
	break;
    
    case GLUT_KEY_PAGE_DOWN: // tilt down
	z += 0.2f;
	lookupdown += 1.0f;
	break;*/

    case GLUT_KEY_UP: // walk forward (bob head)
	xpos += (float)sin(yrot*piover180) * 0.05f; // piover180 --> a coefficient for conversion from degree to radian.
	cout << yrot << endl;
	cout << sin(yrot * piover180) << endl;
	zpos += (float)cos(yrot*piover180) * 0.05f;	
	if (walkbiasangle >= 359.0f)
	    walkbiasangle = 0.0f;	
	else 
	    walkbiasangle+= 10;
	walkbias = (float)sin(walkbiasangle * piover180)/20.0f;
	break;

    case GLUT_KEY_DOWN: // walk back (bob head)
	xpos -= (float)sin(yrot*piover180) * 0.05f;
	zpos -= (float)cos(yrot*piover180) * 0.05f;	
	if (walkbiasangle <= 1.0f)
	    walkbiasangle = 359.0f;	
	else 
	    walkbiasangle-= 10;
	walkbias = (float)sin(walkbiasangle * piover180)/20.0f;
	break;

    case GLUT_KEY_LEFT: // look left
//	yrot += 1.5f;
	xpos += (float)cos(yrot*piover180) * 0.05f; // piover180 --> a coefficient for conversion from degree to radian.
	zpos -= (float)sin(yrot*piover180) * 0.05f;	
	if (walkbiasangle >= 359.0f)
	    walkbiasangle = 0.0f;	
	else 
	    walkbiasangle+= 10;
	walkbias = (float)sin(walkbiasangle * piover180)/20.0f;
	break;

    case GLUT_KEY_RIGHT: // look right
//	yrot -= 1.5f;
	xpos -= (float)cos(yrot*piover180) * 0.05f; // piover180 --> a coefficient for conversion from degree to radian.
	zpos += (float)sin(yrot*piover180) * 0.05f;	
	if (walkbiasangle >= 359.0f)
	    walkbiasangle = 0.0f;	
	else 
	    walkbiasangle+= 10;
	walkbias = (float)sin(walkbiasangle * piover180)/20.0f;
	break;

    default:
	printf ("Special key %d pressed. No action there yet.\n", key);
	break;
    }	
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

    deltaY /= 10;
    deltaX /= 10;

    switch(g_ButtonDown) {
    case GLUT_LEFT_BUTTON:                  // Rotate in the x-y plane
        g_AngleY -= deltaY;
	yrot += deltaY;
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
}

int main(int argc, char **argv) 
{  
    /* load our world from disk */
    SetupWorld();

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
    glutInitWindowSize(640, 480);  

    /* the window starts at the upper left corner of the screen */
    glutInitWindowPosition(0, 0);  

    /* Open a window */  
    window = glutCreateWindow("Computer Graphics Project 3");  

    /* Register the function to do all our OpenGL drawing. */
    glutDisplayFunc(&DrawGLScene);  

    /* Go fullscreen.  This is as soon as possible. */
    glutFullScreen();

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

float *Viewing::ComputeViewMatrix()
{
	Point f;
	f.x = g_view.centerPos ->x - g_view.cameraPos ->x;
	f.y = g_view.centerPos ->y - g_view.cameraPos ->y;
	f.z = g_view.centerPos ->z - g_view.cameraPos ->z;
	f.NormalizeVector();
	Point up;
	up.x = g_view.UP ->x;
	up.y = g_view.UP ->y;
	up.z = g_view.UP ->z;
	up.NormalizeVector();

	Point s;
	Point::CrossProduct(f, up, s);
	s.NormalizeVector();
	Point u;
	Point::CrossProduct(s, f, u);
	u.NormalizeVector();
 
	float * M = new float[16];
	M[0] = s.x; M[1] = s.y; M[2] = s.z; M[3] = 0;
	M[4] = u.x; M[5] = u.y; M[6] = u.z; M[7] = 0;
	M[8] = -f.x; M[9] = -f.y; M[10] = -f.z; M[11] = 0;
	M[12] = 0; M[13] = 0; M[14] = 0; M[15] = 1;

/*	M[0] = s.x; M[4] = s.y; M[8] = s.z; M[12] = 0;
	M[1] = u.x; M[5] = u.y; M[9] = u.z; M[13] = 0;
	M[2] = -f.x; M[6] = -f.y; M[10] = -f.z; M[14] = 0;
	M[3] = 0; M[7] = 0; M[11] = 0; M[15] = 1;*/
	float elem[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, elem);
	for(int i = 0 ; i < 4 ; i++)
	{
		for(int j = 0 ; j < 4 ; j++)
		{
			cout << elem[i * 4 + j] << " ";
		}
		cout << endl;
	}

	return	M;
}


void Point::CrossProduct(const Point p1, const Point p2, Point &p3)
{
	p3.x = p1.y * p2.z - p1.z * p2.y;
	p3.y = p1.z * p2.x - p1.x * p2.z;
	p3.z = p1.x * p2.y - p1.y * p2.x;
}







