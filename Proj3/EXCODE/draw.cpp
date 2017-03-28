/* your drawing code goes in here */

/* by including draw.h, we ensure that the exported prototypes
   match the function definitions */
#include "draw.h"

// Apple's annoying non-standard GL include location
#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <stdlib.h>
#include "example.cpp"
GLfloat verts[] = {50,0,0,  0,50,0,  -50,0,0};
GLfloat colors[] = {1,0,0,  0,1,0,  0,0,1};
GLubyte indices[] = {0, 1, 2};

/*
 * this is called every time the screen needs to be redrawn 
 */

/* called on any keypress
 *
 * We don't use x and y, but they're the mouse location when the key
 * was pressed.
 */
extern "C" void key(unsigned char k, int x, int y)
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
