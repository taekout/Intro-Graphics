/* set up and maintain view as window sizes change */

extern int winWidth, winHeight;
extern float viewWidth;

/* this is called when window is created or resized */
extern "C" void reshape(int width, int height);
