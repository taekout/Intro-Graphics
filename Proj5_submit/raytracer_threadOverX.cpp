#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>

//#define NUM_THREADS	10
int	NUM_THREADS;

using namespace std;
using std::vector;

#define		SPACE_END	100000.0f

#define		FILE_NAME	"balls1.nff"
#define		TRUE	1
#define		FALSE	0
#define		HRESULT		int

#define		CLIPPING(x)	x=(x>255.0)?255.0:x;x=(x<0.0)?0.0:x;
#define		ABS(x)		(x<0.0)?-x:x
#define		SMALLER(x,y)	(x<y)?x:y

/*

unsigned char ***pixels;
pixels = new unsigned char **[50];
for(int i = 0 ; i < 50 ; i ++)
{
	pixels[i] = new unsigned char *[6];
	for(int k = 0 ; k < 6 ; k++)
	{
		pixels[i][k] = new unsigned char[3];
	}
}

*/

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

typedef struct ThreadData
{
	int i;
	int j;
	COLOR3 color;
} ThreadData;

void MatrixMult3X3(float M1[][3], float M2[][3], float resultM[][3]);
void DotProduct(float a[3], float b[3], float *c);
void CrossProduct(float M1[3], float M2[3], float resultM[3]);
void MatrixSubtract1X3(float M1[3], float M2[3], float resultM[3]);
void MatrixSum1X3(float M1[3], float M2[3], float resultM[3]);
void NormalizeVector(float M[3], float resultM[3]);

class Raytracer
{
public:
	// Member Variables
	ifstream	m_file;
	ofstream	m_outFile;
	string		m_output_image_filename;
//	Parser		m_parser;
	VIEW		m_view;
	COLOR3		m_background;
	SHADING		m_shading;

	//Viewing variables
	vector<vector<vector<unsigned char> > >	pixels;

	//Sphere		*m_pSphere;
	list<SPHERE>	m_lSphere;
	int		m_nSphere;
	//Tetra
	list<TETRA>	m_lTetra;
	int		m_nTetra;

	// Member Functions
	Raytracer(string nff_filename, string image_filename);
	~Raytracer();

	// Ray parser Reader Functions
	void Read();
	int OpenFile(string filename);

	// Parsing
	void DoComment();
	void DoViewPoint();
	void DoLightSource();
	void DoBackground();
	void DoFill();
	void DoCone();
	void DoSphere();
	void DoPoly();

	HRESULT ImageAllocation(int nCol, int nRow);	//nCol == nWidth, nRow == nHeight
	void *RayIntersectWithSpheres(void *p);
	HRESULT	Raytracing();
	float GetBarycentricCoordinate(float ray_origin[3], float ray_direction[3]);
};

class Raytracer *g_ray;
static void * ThreadFunc(void *p)
{
	g_ray ->RayIntersectWithSpheres(p);
	pthread_exit((void*) p);
}

int main(int argc, char *argv[])
{
	if(argc == 1)
	{
		cout << "Usage > Program Argument1<Filename>" << endl;
		exit(-1);
	}
	if(argc != 3)
	{
		cout << "There must be 2 parameters" << endl <<
		"Usage>Raytracer.out <NFF file name> <PPM file name>" << endl;
		exit(-2);
	}  
	cout << "How many threads do you want to use? : <Must be an integer #> " << endl;
	cin >> NUM_THREADS;

	// Object Oriented Programming style - Invoke Ray -> Load File.
	string arg1 = argv[1];
	string arg2 = argv[2];
	Raytracer ray(arg1, arg2);
	g_ray = &ray;
	ray.Read();
	ray.Raytracing();


	return	0;
}

Raytracer::Raytracer(string nff_filename, string image_filename)
{
	int result;
	result = OpenFile(nff_filename);
	if(result == FALSE)
	{
		cout << "File does not exist" << endl;
		exit(-1);
	}
	m_nSphere = 0;
	m_nTetra = 0;
	m_output_image_filename = image_filename;
}

Raytracer::~Raytracer()
{
}

void Raytracer::Read()
{
	string line;
	if(!m_file.is_open())
	{
		cout << "File not opened but parser.read() was invoked." << endl;
		exit(-1);
	}
	while(! m_file.eof())
	{
		string sCommand;
	//	getline(m_parser.m_file, line);
	//	cout << line << endl;
		m_file >> sCommand;
	//	cout << sCommand << endl;
		cout << sCommand << endl;
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
        		        break;
        		case 'b':            /* background color */
				DoBackground();
        		        break;
        		case 'f':            /* fill material */
				DoFill();
        		        break;
        		case 'c':            /* cylinder or cone */
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
	// Check data
	/*list<SPHERE>::iterator it;
	for(it = m_lSphere.begin() ; it != m_lSphere.end() ; ++it)
	{
//		cout << m_lSphere.back().x << " " << m_lSphere.back().y << " " << m_lSphere.back().z << " " <<
	//					m_lSphere.back().radius <<  endl;
		cout << (*it).x << " " << (*it).y << " " << (*it).z << " " << (*it).radius << endl;
	}*/
}

int Raytracer::OpenFile(string filename)
{
	m_file.open(filename.c_str());
	if(m_file.is_open() == false)
	{
		return	FALSE;
		cout << "File does not exist" << endl;
		exit(-1);
	}
	return	TRUE;
}

void Raytracer::DoComment()
{
	char sComment[500];
	m_file.getline(sComment, 500);
}

void Raytracer::DoViewPoint()
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
	//printf("%f -- ", m_view.camera.x);
	m_file >> m_view.camera.y;
	//cout << m_view.camera.y << endl;
	//printf("%f -- ", m_view.camera.y);
	m_file >> m_view.camera.z;
	//cout << m_view.camera.z << endl;
	//printf("%f -- \n", m_view.camera.z);

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
	cout << "D NEAR " << m_view.dNear << endl;

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

void Raytracer::DoLightSource()
{
}

void Raytracer::DoBackground()
{
	m_file >> m_background.R;
	m_file >> m_background.G;
	m_file >> m_background.B;
}

void Raytracer::DoFill()
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

void Raytracer::DoCone()
{
}

void Raytracer::DoSphere()
{
	// # of spheres increments
	m_nSphere ++;
	// Now allocate memory for a sphere. But before that, parse it to extract data first.
	SPHERE sp;
	m_file >> sp.x;
	m_file >> sp.y;
	m_file >> sp.z;
	m_file >> sp.radius;

	m_lSphere.push_front(sp);
}

void Raytracer::DoPoly()
{
	int nVertices;
	m_file >> nVertices;
	if(nVertices < 3)
	{
		cout << "Polygons have more than 3 vertices." << endl;
		exit(-6);
	}
	if(nVertices >= 3)
	{
		int nTriangles = nVertices - 3 + 1;
		float vertices[nVertices][3];
		for(int i = 0 ; i < nVertices ; i++)
		{
			m_file >> vertices[i][0];
			m_file >> vertices[i][1];
			m_file >> vertices[i][2];
			cout <<" i = " << i  << "vertex z = "<< vertices[i][2] << endl;
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

HRESULT Raytracer::ImageAllocation(int nCol, int nRow)	//nCol == nWidth, nRow == nHeight
{
	pixels.resize(nRow);
	for(int j = 0 ; j < nRow ; j++)
	{
		pixels[j].resize(nCol);
		for(int i = 0 ; i < nCol ; i++)
		{
			pixels[j][i].resize(3);
		}
	}

	return	1;
}

void *Raytracer::RayIntersectWithSpheres(void *p)
{
	ThreadData *thData = (ThreadData *) p;
	float i = thData ->i; float j = thData ->j;
	//RayIntersect
	//1st compute u,v coordinate
	float Us;
	float Vs;
	//	I, J come from arguments
	// FYI, right - left == resolution.x;
	// u - b == resolution.y; Theoretically, it should not be resolution.y but it does not matter because
	// we are not trying to resize the coordinate. All we need is the equation for the ray, so it is easier
	// just to use the coordinate and resolution.
	// x == half width-length of view plane.
	float x = m_view.dNear * tan(float(m_view.angle) / 180 * 3.141592 / 2); // degree ==> radian!!!
	float r = x;//((float)(m_view.resolution.x - 1)) / 2;
	float l = -x;//-r;
	float t = x;//((float)(m_view.resolution.y - 1)) / 2;
	float b = -x;//-u;
	Us = l + (r - l) * ((float)i + 0.5) / (m_view.resolution.x /*- 1*/);
	Vs = b + (t - b) * ((float)j + 0.5) / (m_view.resolution.y /*- 1*/);
	// I know the intersecting viewplane point (Us,Uv).
	// I need to now know the camera up , camera right, camera view vector.

	// Camera view vector ... get it! - w vector
	// w = lookAt - camera; Dont forget to normalize w.
	float w[3];
	float fLookAt[3], fCameraPos[3];
	fLookAt[0] = m_view.lookAt.x;
	fLookAt[1] = m_view.lookAt.y;
	fLookAt[2] = m_view.lookAt.z;
	fCameraPos[0] = m_view.camera.x;
	fCameraPos[1] = m_view.camera.y;
	fCameraPos[2] = m_view.camera.z;
	MatrixSubtract1X3(fLookAt, fCameraPos, w);
	NormalizeVector(w, w);
	// Camera up vector.. get it! - camera_up[3]
	// To get it, I set an arbitrary vector and then get the other 2 vectors from it.
	float arbitrary_vector[3]= {m_view.cameraUp.x, m_view.cameraUp.y, m_view.cameraUp.z};
	float camera_right[3];
	CrossProduct(w, arbitrary_vector, camera_right);
	NormalizeVector(camera_right, camera_right);
	float camera_up[3];
	CrossProduct(camera_right, w, camera_up);
	NormalizeVector(camera_up, camera_up);

	// Now we know camera_up, w(viewing direction), camera_right)
	// s = e + Us*U + Vs*V - d*W (s = origin(e) + direction(Us*U + Vs*V - d*w))
	float s[3];
	float ray_direction[3] = {0.0f, 0.0f, 0.0f};// = {camera_right, camera_up, w};
	camera_right[0] *= Us; camera_right[1] *= Us; camera_right[2] *= Us;
	// Image plane starts from top, which is (0,0) to bottom(512,512)
	camera_up[0] *= -Vs; camera_up[1] *= -Vs; camera_up[2] *= -Vs;
	MatrixSum1X3(camera_right, camera_up, ray_direction);
	w[0] *= m_view.dNear; w[1] *= m_view.dNear; w[2] *= m_view.dNear;
	MatrixSum1X3(ray_direction, w, ray_direction);
	NormalizeVector(camera_right, camera_right);
	NormalizeVector(camera_up, camera_up);
	NormalizeVector(w, w);
	float ray_origin[3];
	ray_origin[0] = m_view.camera.x; ray_origin[1] = m_view.camera.y; ray_origin[2] = m_view.camera.z;
	MatrixSum1X3(ray_direction, ray_origin, s);
	
	//  Now we have Us, Vs, S, so just compute t right away.
	// Compute parameter t for all the spheres in LIST.
	// and then get the closest t. We are going to use the closest t to display colors.
	float d[3];
	float e[3];
	e[0] = m_view.camera.x; e[1] = m_view.camera.y; e[2] = m_view.camera.z;
	// d = s - e;
	d[0] = s[0] - e[0]; d[1] = s[1] - e[1]; d[2] = s[2] - e[2];
	float d_dot_d;
	DotProduct(d, d, &d_dot_d);
	float minus_d[3];
	minus_d[0] = -d[0]; minus_d[1] = -d[1]; minus_d[2] = -d[2];
	// e_c varies for each sphere. so I must compute e_c for all the spheres.
	// Radius also varies for each sphere.
	// Conclusion ! > I must compute R, e_c for each sphere.
	float e_c[3];
	float radius;
	float c[3];
	//list<SPHERE> it;
	float closest_t = SPACE_END;
	float closest_c[3] = {SPACE_END, SPACE_END, SPACE_END};
	float closest_R = SPACE_END;
	for(list<SPHERE>::iterator it = m_lSphere.begin() ; it != m_lSphere.end() ; it++)
	{
		// c = center of each sphere
		c[0] = (*it).x;	c[1] = (*it).y; c[2] = (*it).z;
		radius = (*it).radius;
		e_c[0] = e[0] - c[0]; e_c[1] = e[1] - c[1]; e_c[2] = e[2] - c[2];
		float determinant = 0.0f;
		float tmp1;
		DotProduct(d, e_c, &tmp1);
		tmp1 = tmp1 * tmp1;
		float tmp2;
		DotProduct(e_c, e_c, &tmp2);
		tmp2 -= (radius * radius);
		float tmp3;
		tmp3 = d_dot_d * tmp2;
		determinant = tmp1 - tmp3;
		// Just for debugging. I am wary of floating point precision issues because it might not be 0 but 0.0000001
//		if(0.0000001 > determinant && determinant > -0.0000001)
//		{
//			cout << "There might be precision issue. Check line 463" << endl;
//		}
		if(t < 0)		// t < 0 ==> The sphere is behind the camera ==> Does not Intersect!
			continue;
		if(determinant > 0)	// 2 roots
		{
			float t1, t2;
			float term1;
			DotProduct(minus_d, e_c, &term1);
			t1 = term1 + determinant;
			t1 /= d_dot_d;
			t2 = term1 - determinant;
			t2 /= d_dot_d;
			if(closest_t > t1)
			{	// Whenever we find the t, compute N and save it in the global variable.
				closest_t = t1;
				closest_c[0] = c[0];closest_c[1] = c[1];closest_c[2] = c[2];
				closest_R = radius;
			}
			//closest_t = SMALLER(closest_t, t1);
			if(closest_t > t2)
			{
				closest_t = t2;
				closest_c[0] = c[0];closest_c[1] = c[1];closest_c[2] = c[2];
				closest_R = radius;
			}
			//closest_t = SMALLER(closest_t, t2);
		}
		else if(determinant == 0)	// only 1 root
		{
			float t;
			DotProduct(minus_d, e_c, &t);
			t /= d_dot_d;
			if(closest_t > t)
			{
				closest_t = t;
				closest_c[0] = c[0];closest_c[1] = c[1];closest_c[2] = c[2];
				closest_R = radius;
			}
			//closest_t = SMALLER(closest_t, t);
		}
		else				// no roots -- Not intersects!
		{
			;
		}

	}
	float t_fromTetra = GetBarycentricCoordinate(ray_origin, ray_direction);
	if(t_fromTetra < closest_t)
		closest_t = t_fromTetra;
	// At this point, we know the closest t. If t == 100000.0(MACRO:SPACE_END), the ray does NOT intersect with any spheres.
	if(closest_t == SPACE_END)
	{
		// Never intersect with any spheres that we have.
		// In this case, just return color = {-1, -1, -1};
		//return color;
		thData ->color.R = -1;
		thData ->color.G = -1;
		thData ->color.B = -1;
		pthread_exit((void *)thData);
	}

	// It is time to compute P(t) = e + t(s-e)
	float P[3];
	P[0] = e[0] + closest_t * (s[0] - e[0]);
	P[1] = e[1] + closest_t * (s[1] - e[1]);
	P[2] = e[2] + closest_t * (s[2] - e[2]);
	// Now we know P, Center and radius of the closest sphere, which we need to compute N.
	// so compute N.
	float N[3];
	N[0] = (P[0] - closest_c[0]) / closest_R;
	N[1] = (P[1] - closest_c[1]) / closest_R;
	N[2] = (P[2] - closest_c[2]) / closest_R;

	if(closest_t != SPACE_END)
	{
		thData ->color.R = m_shading.color.R;
		thData ->color.G = m_shading.color.G;
		thData ->color.B = m_shading.color.B;
		//return	color;
		pthread_exit((void *) thData);
	}
	//return	color;
	thData ->color.R = -2;
	thData ->color.G = -2;
	thData ->color.B = -2;
	pthread_exit((void *) thData);
}

HRESULT	Raytracer::Raytracing()
{
//	Here, Use thread or whatever.. Just traverse with i, j.
//	i = [l, r], r = [t,b]
	// Now save the data into the vectors(== array)
	pthread_t thread[NUM_THREADS];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	int rc;
	int resY = m_view.resolution.y , resX = m_view.resolution.x;
	ImageAllocation(m_view.resolution.x, m_view.resolution.y);
	ThreadData thData[NUM_THREADS];
	struct timeval timest1, timest2;
	gettimeofday(&timest1, NULL);
	for(int j = 0 ; j < resY ; j++)
	{
//		cout << "Debugging code : " << j << endl;
		int nActiveThread = 0;
		for(int i = 0 ; i < resX ; i += nActiveThread)
		{
//cout << "Active Thread # : " << nActiveThread << endl;
			COLOR3 color;
			nActiveThread = min(NUM_THREADS, resX - i);
			for(int k = 0 ; k < nActiveThread ; k++)
			{
				thData[k].i = i + k; thData[k].j = j;
				rc = pthread_create(&thread[k], &attr, ThreadFunc, (void *)&thData[k]);
				if(rc != 0)
				{
					cout << "Thread Creation failed. Error code #" << rc << endl;
					exit(-7);
				}
			}
			// Thread Join
			int *nStatus;		
			for(int k = 0 ; k < nActiveThread ; k++)
			{
				rc = pthread_join(thread[k], (void **)&nStatus);
				if(rc != 0)
				{
					cout << "Thread Join failed. Error Code # " << rc << endl;
				}
			}

			//color = RayIntersectWithSpheres(i, j);
			for(int k = 0 ; k < nActiveThread ; k++)
			{
				COLOR3 color;
				color.R = thData[k].color.R;
				color.G = thData[k].color.G;
				color.B = thData[k].color.B;
//cout << "R : " << color.R << " G : " << color.G << " B : " << color.B << endl;
				if(color.R == -1 && color.G == -1 && color.B == -1)
				{
					// Never intersect
					// Later, in this part, color should be the background color.
					color.R = m_background.R;
					color.G = m_background.G;
					color.B = m_background.B;
//cout << "BACKGROUND COLOR SET" << endl;
				}
				color.R = (color.R * 255);
				color.G = (color.G * 255);
				color.B = (color.B * 255);
				CLIPPING(color.R);
				CLIPPING(color.G);
				CLIPPING(color.B);
				unsigned char R = (unsigned char)color.R;
				unsigned char G = (unsigned char)color.G;
				unsigned char B = (unsigned char)color.B;
//cout << "Actual color: R-" << R << " G-" << G << " B-" << B << endl;
				pixels[j][i+k][0] = R;
				pixels[j][i+k][1] = G;
				pixels[j][i+k][2] = B;
			}
		}
	}
	gettimeofday(&timest2, NULL);
	long mtime, seconds, useconds;
	seconds = timest2.tv_sec - timest1.tv_sec;
	useconds = timest2.tv_usec - timest1.tv_usec;
	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
	cout << "Elapsed time : " << mtime << endl;
	m_outFile.open(m_output_image_filename.c_str());
	unsigned char buffer[resX * 3];
	m_outFile << "P6\n" << m_view.resolution.x << " " << m_view.resolution.y << "\n255\n";
	for(int j = 0 ; j < resY ; j++)
	{
		for(int i = 0 ; i < resX ; i++)
		{
			buffer[i * 3/*com size*/ /*+ 0*/] = pixels[j][i][0];
			buffer[i * 3/*com size*/ + 1] = pixels[j][i][1];
			buffer[i * 3/*com size*/ + 2] = pixels[j][i][2];
		}
		m_outFile.write((char *)buffer, resX * 3 * sizeof(unsigned char));
	}
}


void MatrixMult3X3(float M1[][3], float M2[][3], float resultM[][3])
{
	int i = 0, j = 0 , k = 0;
	for(i = 0 ; i < 2 ; i++)
	{
		for(j = 0 ; j < 3 ; j++)
		{
			resultM[i][j] = 0.0f;
			for(k = 0 ; k < 3 ; k++)
			{
				resultM[i][j] = M1[i][k] * M2[k][j];
			}
		}
	}
}

void DotProduct(float a[3], float b[3], float *result)
{
	*result = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void MatrixSubtract1X3(float M1[3], float M2[3], float resultM[3])
{
	resultM[0] = M1[0] - M2[0];
	resultM[1] = M1[1] - M2[1];
	resultM[2] = M1[2] - M2[2];
}

void MatrixSum1X3(float M1[3], float M2[3], float resultM[3])
{
	resultM[0] = M1[0] + M2[0];
	resultM[1] = M1[1] + M2[1];
	resultM[2] = M1[2] + M2[2];
}


void NormalizeVector(float M[3], float resultM[3])
{
	float length = M[0] * M[0] + M[1] * M[1] + M[2] * M[2];
	length = sqrt(length);
	resultM[0] = M[0] / length;
	resultM[1] = M[1] / length;
	resultM[2] = M[2] / length;
}


void CrossProduct(float M1[3], float M2[3], float resultM[3])
{
	resultM[0] = M1[1] * M2[2] - M1[2] * M2[1];
	resultM[1] = M1[2] * M2[0] - M1[0] * M2[2];
	resultM[2] = M1[0] * M2[1] - M1[1] * M2[0];
}


float Raytracer::GetBarycentricCoordinate(float ray_origin[3], float ray_direction[3])
{
	float closest_t = SPACE_END;
	for(list<TETRA>::iterator it = m_lTetra.begin() ; it != m_lTetra.end() ; it++)
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
			cout << "The triangle is degenerate or parallel to the ray." << endl;
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
		}
	}
	return	closest_t;
}

/*
#include <pthread.h>

#define NUM_THREADS     5

void *PrintHello(void *threadid)
{
   long tid;
   tid = (long)threadid;
   printf("Hello World! It's me, thread #%ld!\n", tid);
   pthread_exit(NULL);
}

int main (int argc, char *argv[])
{
// MUST COMPILE with g++     ============+> with -lpthread or -pthread option!!!!!
   pthread_t threads[NUM_THREADS];
   int rc;
   long t;
   for(t=0; t<NUM_THREADS; t++){
      printf("In main: creating thread %ld\n", t);
      rc = pthread_create(&threads[t], NULL, PrintHello, (void *)t);
      if (rc){
         printf("ERROR; return code from pthread_create() is %d\n", rc);
         exit(-1);
      }
   }
   pthread_exit(NULL);
}
*/
