#include <gl/glut.h>
#include <math.h>

float g_aspect;
float g_time = 0;
float g_dt = 0.05;

int g_w;
int g_h;

float g_l[3] = {0,5,5};			// 조명의 위치를 설정
float g_v[3] = {0,0,5};			// 카메라의 위치를 설정

void ComputeColor(float l[], float v[], float n[],
				  float c[], float ka, float kd, float ks)
{
	// Phong Illumination Implementation



}

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

	// 점들의 좌표 계산
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

	// 맨 마지막 점
	v[vcount][0] = -r;
	v[vcount][1] = 0;
	v[vcount][2] = 0;
	vcount++;

	glBegin(GL_POINTS);
	for(i=0; i<num_vertex; i++)
		glVertex3fv(v[i]);
	glEnd();


	// 뚜껑
	int fcount = 0;
	for(int i=0; i<div; i++)
	{
		f[fcount][0] = 0;
		f[fcount][1] = i+1;
		if(i<div-1)
			f[fcount][2] = i+2;
		else		// 마지막 점
			f[fcount][2] = 1;
		fcount++;
	}
	
	// 몸체
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
		float p[3]={};									// 면의 대표점
		float l[3]={};									// 빛의 방향 (조명의 방향)
		float vc[3]={};									// 카메라 방향
		float n[3]={};									// 법선벡터
		float ka = 0.2;									// 보통 이정도 준다. (합하면 1이 되도록)
		float kd = 0.5;
		float ks = 0.3;
		float c[3]={};									// 색깔은 RGB니까 3개

		for(int k=0; k<3; k++)							// k는 면에서의 점1, 점2, 점3
			for(int j=0; j<3; j++)						// j는 x,y,z의 인덱스
			{
				p[j] += v[f[i][k]][j];					// p의 x,y,z는 어떤 면의 3점(f의 i번째 면의 k번째 점)의 합
			}

		for(int j=0; j<3; j++)
			p[j]/= 3.0f;


		for(int j=0; j<3; j++)
		{
			l[j] = g_l[j]-p[j];							// 이점에서 바라본 빛의 위치는
			vc[j] = g_v[j]-p[j];						// 이점에서 바라본 카메라의 위치는
			n[j] = p[j];								// 이점에서의 노멀방향은 (구이기때문에 p와 동일함)
		}
		float leng1, leng2, leng3;						// 우리는 필요한게 unit vector(단위벡터)다 (방향벡터이기 때문에)
		leng1 = sqrt(l[0]*l[0]+l[1]*l[1]+l[2]*l[2]);
		leng2 = sqrt(vc[0]*vc[0]+vc[1]*vc[1]+vc[2]*vc[2]);
		leng3 = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
		for(int j=0; j<3; j++)
		{
			l[j] /= leng1;
			vc[j] /= leng2;
			n[j] /= leng3;
		}
		
		ComputeColor(l, vc, n, c, ka, kd, ks);			// 이 함수를 부르면 칼라값이 나오도록

		glColor3fv(c);									// 이 색깔을 주면 된다.
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
	glEnable(GL_DEPTH_TEST);
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
	float y = 0;
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