#include <cstdlib>

#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#elif defined(_WIN32)
#include "glut.h"
#else
#include <GL/glut.h>
#endif

// A set of vertices for a cube.
static GLfloat g_Vertices[] = {
     0.500000f,  0.500000f,  0.500000f,
    -0.500000f,  0.500000f,  0.500000f,
     0.500000f, -0.500000f,  0.500000f,
    -0.500000f, -0.500000f,  0.500000f,
     0.500000f,  0.500000f, -0.500000f,
    -0.500000f,  0.500000f, -0.500000f,
     0.500000f, -0.500000f, -0.500000f,
    -0.500000f, -0.500000f, -0.500000f
};

// A set of normals for the vertices.
static GLfloat g_Normals[] = {
     0.577350f,  0.577350f,  0.577350f,
    -0.333333f,  0.666667f,  0.666667f,
     0.666667f, -0.333333f,  0.666667f,
    -0.666667f, -0.666667f,  0.333333f,
     0.666667f,  0.666667f, -0.333333f,
    -0.666667f,  0.333333f, -0.666667f,
     0.333333f, -0.666667f, -0.666667f,
    -0.577350f, -0.577350f, -0.577350f
};

// A set of vertex colors for the vertices.
static GLfloat g_Colors[] = {
    0.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 1.0f,
    0.0f, 1.0f, 0.0f,
    0.0f, 1.0f, 1.0f,
    1.0f, 0.0f, 0.0f,
    1.0f, 0.0f, 1.0f,
    1.0f, 1.0f, 0.0f,
    1.0f, 1.0f, 1.0f
};

// A set of indices for 4 triangles to form a cube.
static GLuint g_Indices[] = {
    0, 1, 2,
    3, 2, 1,
    0, 2, 4,
    6, 4, 2,
    0, 4, 1,
    5, 1, 4,
    7, 5, 6,
    4, 6, 5,
    7, 6, 3,
    2, 3, 6,
    7, 3, 5,
    1, 5, 3
};

void drawCube() {
    // Specify the vertex, color, and normal data.
    glVertexPointer(3, GL_FLOAT, 0, g_Vertices);
    glEnableClientState(GL_VERTEX_ARRAY);

    glColorPointer(3, GL_FLOAT, 0, g_Colors);
    glEnableClientState(GL_COLOR_ARRAY);

    glNormalPointer(GL_FLOAT, 0, g_Normals);
    glEnableClientState(GL_NORMAL_ARRAY);


    // Draw the specified triangles.
    int triangles = sizeof(g_Indices) / sizeof(g_Indices[0]);
    glDrawElements(GL_TRIANGLES, triangles, GL_UNSIGNED_INT, g_Indices);


    // Reset what we turned on.
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
}
