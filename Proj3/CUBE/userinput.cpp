#include <cstdlib>

#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#elif defined(_WIN32)
#include "glut.h"
#else
#include <GL/glut.h>
#endif

#include "userinput.h"
#include <iostream>

using namespace std;

bool g_ShiftDown;
int g_ButtonDown, g_LastX = -1, g_LastY = -1;
double g_AngleX = 0.0, g_AngleY = 0.0, g_AngleZ = 0.0;
double g_PositionX = 0.0f, g_PositionY = 0.0, g_PositionZ = -4.0;

// Handle any keyboard events.
// key contains the ASCII code of the key pressed.
void keyboard(unsigned char key, int x, int y) {
    // Exit on escape.
    if (key == 27) exit(0);
	if(key == 32) cout << "SPACE " << endl;//printf("SPACE");
}

// Handle any mouse button events.
void mouse(int button, int state, int x, int y) {
    // Record the mouse position when a button is pressed.
    g_LastX = x;
    g_LastY = y;
    g_ButtonDown = button;
    g_ShiftDown = (glutGetModifiers() == GLUT_ACTIVE_SHIFT);
}

// Handle any mouse movement events.
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
}
