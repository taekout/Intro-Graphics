#include <gl/glut.h>
#include <math.h>

float g_aspect;
float g_time = 0;
float g_dt = 0.05;

int g_w;
int g_h;


void MySphere(float r, int div)
{
	if(div<2) return;

	float ** v;
	int num_vertex = (div-1)*div + 2;
	v = new float * [num_vertex];
	int i,j;
	for(i=0; i<num_vertex; i++)
		v[i] = new float[3];

	int ** f;
	int num_face = div*2 + div*2*(div-1);
	f = new int * [num_face];
	for(i=0; i<num_face; i++)
		f[i] = new int[3];

	
	int vcount = 0;
	v[vcount][0] = r;
	v[vcount][1] = 0;
	v[vcount][2] = 0;
	vcount++;

	float dt1 = 180.0f/div* 3.141592/180.0f;
	float dt2 = 360.0f/div* 3.141592/180.0f;

	for(int i=1; i<div; i++)
	{
		float th1 = dt1*i;
		float x = r*cos(th1);
		float r1 = r*sin(th1);

		for(int j=0; j<div; j++)
		{
			float th2 = dt2*j;
			float y = r1*cos(th2);
			float z = r1*sin(th2);
			v[vcount][0] = x;
			v[vcount][1] = y;
			v[vcount][2] = z;
			vcount ++;
		}
	}

	v[vcount][0] = -r;
	v[vcount][1] = 0;
	v[vcount][2] = 0;
	vcount++;

	glBegin(GL_POINTS);
	for(i=0; i<num_vertex; i++)
		glVertex3fv(v[i]);
	glEnd();


	int fcount = 0;
	for(int i=0; i<div; i++)
	{
		f[fcount][0] = 0;
		f[fcount][1] = i+1;
		if(i<div-1)
			f[fcount][2] = i+2;
		else
			f[fcount][2] = 1;
		fcount++;
	}
	
	for(int i=0; i<div-2; i++)
	{
		int v_st = div*i + 1;
		
		for(int j=0; j<div; j++)
		{
			if(j<div-1)
			{
				f[fcount][0] = v_st+j;
				f[fcount][1] = v_st+div+j;
				f[fcount][2] = v_st+div+j+1;
				fcount ++;
				f[fcount][0] = v_st+div+j+1;
				f[fcount][1] = v_st+j+1;
				f[fcount][2] = v_st+j;
				fcount ++;
			}
			else
			{
				f[fcount][0] = v_st+j;
				f[fcount][1] = v_st+div+j;
				f[fcount][2] = v_st+div;
				fcount ++;
				f[fcount][0] = v_st+div;
				f[fcount][1] = v_st;
				f[fcount][2] = v_st+j;
				fcount ++;
			}
		}
	}

	int v_st = div*(div-2)+1;
	for(int i=0; i<div; i++)
	{
		f[fcount][0] = v_st+i;
		f[fcount][1] = div*(div-1)+1;
		if(i<div-1)
			f[fcount][2] = v_st+i+1;
		else
			f[fcount][2] = v_st;
		fcount++;
	}

	glBegin(GL_TRIANGLES);
	for(i=0; i<fcount; i++)
	{
		glVertex3fv(v[f[i][0]]);
		glVertex3fv(v[f[i][1]]);
		glVertex3fv(v[f[i][2]]);
	}
	glEnd();


	for(i=0; i<num_vertex; i++)
		delete [] v[i];
	delete [] v;

	for(i=0; i<num_face; i++)
		delete [] f[i];
	delete [] f;


}


void MyInit()
{
	glClearColor(0,0,0,0);
	
//	glEnable(GL_LIGHTING);
//	glEnable(GL_LIGHT0);
//	glEnable(GL_COLOR_MATERIAL);
//	glEnable(GL_DEPTH_TEST);
}

void MyIdle()
{
	g_time += g_dt;
	glutPostRedisplay();
}




void MyReshape(int w, int h)
{
	g_w = w;
	g_h = h;
	
	glViewport(0,0,w,h);
	g_aspect = float(w)/float(h);
}
void MyDraw()
{
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	
	////////////////////////////////////////////////
	//	parallel projection

	//	glOrtho(-1,1,-1/g_aspect, 1/g_aspect, 1, 10);
	
	////////////////////////////////////////////////
	//	Perspective projection by using glFrustum

	//	float size = 1.0f;
	//	glFrustum(-size, size, -size/g_aspect, size/g_aspect, 1, 10);	

	////////////////////////////////////////////////
	//	Perspective projection by using glPerspective

	float angle = 30.0f;
	gluPerspective(angle, g_aspect, 1, 10);

	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	float rad = g_time;
	float x = sin(rad)*5.0f;
	float z = cos(rad)*5.0f;
	gluLookAt(x,0,z, 0,0,0, 0,1,0);

	glDisable(GL_LIGHTING);
	float leng = 1.4f;
	glBegin(GL_LINES);
		glColor3f(1,0,0);	
		glVertex3f(0,0,0);
		glVertex3f(leng,0,0);
		glColor3f(0,1,0);
		glVertex3f(0,0,0);
		glVertex3f(0,leng,0);
		glColor3f(0,0,1);
		glVertex3f(0,0,0);
		glVertex3f(0,0,leng);
	glEnd();


//	glEnable(GL_LIGHTING);	
	glColor3f(0,0,1);
	MySphere(1, 10);

	glutSwapBuffers();
	glFlush();
}



int main()
{

	
	glutInitDisplayMode(GLUT_DOUBLE);
	
	glutInitWindowSize(400,300);
	glutInitWindowPosition(100,100);
	glutCreateWindow("3D Sphere");

	MyInit();
	glutReshapeFunc(MyReshape);
	glutIdleFunc(MyIdle);
	glutDisplayFunc(MyDraw);

	glutMainLoop();

	return 0;
}