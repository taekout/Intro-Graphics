#ifndef _USERINPUT_H_
#define _USERINPUT_H_

// Keep track of the mouse position for model rotation and translation.
extern bool g_ShiftDown;
extern int g_ButtonDown, g_LastX, g_LastY;
extern double g_AngleX, g_AngleY, g_AngleZ;
extern double g_PositionX, g_PositionY, g_PositionZ;

// Handle any keyboard events.
extern "C" void keyboard(unsigned char key, int x, int y);

// Handle any mouse button events.
extern "C" void mouse(int button, int state, int x, int y);

// Handle any mouse movement events.
extern "C" void motion(int x, int y);

#endif // _USERINPUT_H_
