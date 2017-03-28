#ifndef __DRAW__
#define __DRAW__

/* called for every screen draw */
extern "C" void draw(void);

/* called on any keypress */
extern "C" void key(unsigned char k, int x, int y);

#endif



