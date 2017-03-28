/*
http://ftp.gwdg.de/pub/grafik/vogl/nff/nff_files/
Where I can get some nff files.

I studied and used Blinn Phong from : http://en.wikipedia.org/wiki/Blinn%E2%80%93Phong_shading_model
*/

#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <math.h>
#include <pthread.h>
#include <sys/time.h>
#include <string.h>

#define PARAMETER_MODE
//#define NORMAL_MODE
#define GL_CONSTANT_ATTENUATION  0x1207
#define GL_LINEAR_ATTENUATION    0x1208
#define GL_QUADRATIC_ATTENUATION 0x1209

using namespace std;
using std::vector;
int nSet = 0;
#define		SPACE_END	100000.0f

#define		FILE_NAME	"balls1.nff"
#define		TRUE	1
#define		FALSE	0
#define		HRESULT		int

#define		CLIPPING(x)	x=(x>255.0)?255.0:x;x=(x<0.0)?0.0:x;
#define		ABS(x)		(x<0.0)?-x:x
#define		SMALLER(x,y)	(x<y)?x:y

int g_nActiveThread;


typedef struct COLOR3
{
	float R;
	float G;
	float B;
	COLOR3()
	{
		R = G = B = -1;
	}
	COLOR3(float r, float g, float b)
	{
		R = r ; G = g ; B = b;
	}
	COLOR3 &operator = (const COLOR3 &src)
	{
		this ->R = src.R;
		this ->G = src.G;
		this ->B = src.B;
		return *this;
	}
} COLOR3;

class Point{
public:
	double x;
	double y;
	double z;

public:
	void SanityCheck()
	{
	}
	Point(){x = -1; y = -1; z = -1;}
	Point(double x, double y, double z)
	{
		this ->x = x; this ->y = y; this ->z = z;
	}
	Point(double xyz[3])
	{
		x = xyz[0]; y = xyz[1]; z = xyz[2];
	}
	Point(float xyz[3])
	{
		x = xyz[0]; y = xyz[1]; z = xyz[2];
	}
	Point(const Point &pt)
	{
		x = pt.x ; y = pt.y ; z = pt.z ;
	}
	Point(const COLOR3 &clr)
	{
		x = clr.R; y = clr.G; z = clr.B;
	}
	Point &operator =(const Point & src)
	{
		x = src.x;
		y = src.y;
		z = src.z;
		return *this;
	}
	Point &operator =(const double *src)
	{
		x = src[0];
		y = src[1];
		z = src[2];
		return *this;
	}
	Point &operator =(const float *src)
	{
		x = src[0];
		y = src[1];
		z = src[2];
		return *this;
	}

// ADD * operator!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Point & operator *(double t)
	{
		Point *pt = new Point(this ->x, this ->y, this ->z);
		pt ->x = pt ->x * t;
		pt ->y = pt ->y * t;
		pt ->z = pt ->z * t;
		return *pt;
	}
	Point & operator *(const Point &src)
	{
		Point *pt = new Point(this ->x, this ->y, this ->z);
		pt ->x = pt ->x * src.x;
		pt ->y = pt ->y * src.y;
		pt ->z = pt ->z * src.z;
		return	 *pt;
	}
	Point & operator /(const Point &src)
	{
		Point *pt = new Point(this ->x, this ->y, this ->z);
		pt ->x = pt ->x / src.x;
		pt ->y = pt ->y / src.y;
		pt ->z = pt ->z / src.z;
		return	*pt;
	}
	Point & operator*=(const double t)
	{
		x *= t;
		y *= t;
		z *= t;
		return (*this);
	}
	Point & operator/=(const double t)
	{
		x /= t;
		y /= t;
		z /= t;
		return (*this);
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
		return Point(*this) += src;
	}

	double* GetPoint(void)
	{
//		Point ret(x, y, z);
//		return	ret;
		double  *p = new double[3];
		p[0] = x;
		p[1] = y;
		p[2] = z;
		return p;
	}
	void NormalizeVector()
	{
		double len = GetLengthVector();
		double reciprocalLen = 1.0 / len;
		x *= reciprocalLen;
		y *= reciprocalLen;
		z *= reciprocalLen;
	}
	double GetLengthVector()
	{
		double len = x * x + y * y + z * z;
		len = sqrt(len);
		return len;
	}

	static void CrossProduct(Point &p1, Point &p2, Point &p3)
	{
		double *M1 = p1.GetPoint();
		double *M2 = p2.GetPoint();
		p3.x = M1[1] * M2[2] - M1[2] * M2[1];
		p3.y = M1[2] * M2[0] - M1[0] * M2[2];
		p3.z = M1[0] * M2[1] - M1[1] * M2[0];
	}

	static double DotProduct(const Point &p1, const Point &p2)
	{
		return	p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
	}
	static int SimilarityCheck(Point &a, Point &b)
	{
		double dif = a.x - b.x;
		double dif2 = a.y - b.y;
		double dif3 = a.z - b.z;
		if((dif + dif2 + dif3) / 3.0 < 0.1)
			return	1;
		else
			return	0;
	}
	void Print()
	{
		cout << x << " " << y << " " << z << endl;
	}
	int IsZero()
	{
		if(x == 0 && y == 0 && z == 0)
			return	1;
		else	return	0;
	}
};


typedef struct VIEW
{
	Point	camera;
        Point	lookAt;
        Point	cameraUp;
        float	angle;
        float	dNear;		// In this case, the view plane(projection plane) lies at the near Plane
        struct	RESOLUTION
	{
		int x;
		int y;
	} resolution;
	int nLights;
	struct lightSources
	{
		Point light;
		COLOR3 lightColor;
		float lightRadius;
	} lights[20];
	VIEW()
	{
		nLights = 0;
		for(int i = 0 ; i < 20 ; i ++)
		{
			lights[i].light.x = SPACE_END; lights[i].light.y = SPACE_END; lights[i].light.z = SPACE_END;
			lights[i].lightColor.R = -1; lights[i].lightColor.G = -1; lights[i].lightColor.B = -1;
			lights[i].lightRadius = -1;
		}
	}
} VIEW;

typedef struct SHADING
{
	COLOR3	color;
	float	Kd;
	float	Ks;
	float	PhongPow;
	float	Transmission;
	float	indexRefraction;
	SHADING &operator =(const SHADING &src)
	{
		this ->color.R = src.color.R;
		this ->color.G = src.color.G;
		this ->color.B = src.color.B;
		this ->Kd = src.Kd;
		this ->Ks = src.Ks;
		this ->PhongPow = src.PhongPow;
		this ->Transmission = src.Transmission;
		this ->indexRefraction = src.indexRefraction;
		return *this;
	}
} SHADING;

typedef struct SPHERE
{
	float x;
	float y;
	float z;
	float radius;
	int ID;	// ID starts from [1, 2, 3, ... ]
	SHADING shd;
	SPHERE &operator =(const SPHERE & src)
	{
		this ->x = src.x;
		this ->y = src.y;
		this ->z = src.z;
		this ->radius = src.radius;
		this ->ID = src.ID;
		this ->shd = src.shd;
		return *this;
	} 
} SPHERE;

typedef struct TETRA
{
	double a[3];
	double b[3];
	double c[3];
	double N[3];
	double D[3];
	int ID; // ID starts from [1,2,3, ...]
	SHADING shd;
	TETRA & operator =(const TETRA &src)
	{
		for(int i = 0 ; i < 3 ; i++)
		{
			this ->a[i] = src.a[i];
			this ->b[i] = src.b[i];
			this ->c[i] = src.c[i];
			this ->N[i] = src.N[i];
			this ->D[i] = src.D[i];
		}
		this ->ID = src.ID;
		this ->shd = src.shd;
		return	*this;
	}
} TETRA;

typedef struct TETRA_PP
{
	double a[3];
	double b[3];
	double c[3];
	double a_N[3];
	double b_N[3];
	double c_N[3];
	double D[3];
	int ID;
	SHADING shd;
	TETRA_PP & operator =(const TETRA_PP &src)
	{
		for(int i = 0 ; i < 3 ; i++)
		{
			this ->a[i] = src.a[i];
			this ->b[i] = src.b[i];
			this ->c[i] = src.c[i];
			this ->a_N[i] = src.a_N[i];
			this ->b_N[i] = src.b_N[i];
			this ->c_N[i] = src.c_N[i];
			this ->D[i] = src.D[i];
		}
		this ->ID = src.ID;
		this ->shd = src.shd;
		return	*this;
	}
} TETRA_PP;

void MatrixMult3X3(float M1[][3], float M2[][3], float resultM[][3]);
void DotProduct(float a[3], float b[3], float *c);
void DotProduct(double a[3], double b[3], double *c);
void CrossProduct(float M1[3], float M2[3], float resultM[3]);
void CrossProduct(double M1[3], double M2[3], double resultM[3]);
void MatrixSubtract1X3(float M1[3], float M2[3], float resultM[3]);
void MatrixSubtract1X3(double M1[3], double M2[3], double resultM[3]);
void MatrixSum1X3(float M1[3], float M2[3], float resultM[3]);
void MatrixSum1X3(double M1[3], double M2[3], double resultM[3]);
void NormalizeVector(float M[3], float resultM[3]);
void NormalizeVector(double M[3], double resultM[3]);


typedef struct INTERSECT_SPHERE_INFO
{
	SPHERE sp;
	Point p;
	double t;
	INTERSECT_SPHERE_INFO()
	{
		p = Point(SPACE_END, SPACE_END, SPACE_END);
		t = SPACE_END;
		sp.x = SPACE_END; sp.y = SPACE_END; sp.z = SPACE_END;
		sp.radius = -1; sp.ID = -1;
		sp.shd.color.R = -1; sp.shd.color.G = -1; sp.shd.color.B = -1;
	}
	INTERSECT_SPHERE_INFO & operator =(const INTERSECT_SPHERE_INFO &src)
	{
		this ->sp = src.sp;
		this ->p = src.p;
		this ->t = src.t;
	}
} INTERSECT_SPHERE_INFO;

typedef struct INTERSECT_TRI_INFO
{
	TETRA tet;
	Point p;
	float t;
	INTERSECT_TRI_INFO()
	{
		p = Point(SPACE_END, SPACE_END, SPACE_END);
		t = SPACE_END;
		tet.ID = -1;
		for(int i = 0 ; i < 3 ; i++)
		{
			tet.a[i] = SPACE_END;
			tet.b[i] = SPACE_END;
			tet.c[i] = SPACE_END;
			tet.N[i] = SPACE_END;
			tet.D[i] = SPACE_END;
		}
	}
	INTERSECT_TRI_INFO & operator =(const INTERSECT_TRI_INFO &src)
	{
		this ->tet = src.tet;
		this ->p = src.p;
		this ->t = src.t;
		return	*this;
	}
} INTERSECT_TRI_INFO;

typedef struct INTERSECT_TRI_P_INFO
{
	TETRA_PP tet;
	Point p;
	float t;
	INTERSECT_TRI_P_INFO()
	{
		p = Point(SPACE_END, SPACE_END, SPACE_END);
		t = SPACE_END;
		tet.ID = -1;
		for(int i = 0 ; i < 3 ; i++)
		{
			tet.a[i] = SPACE_END;
			tet.b[i] = SPACE_END;
			tet.c[i] = SPACE_END;
			tet.a_N[i] = SPACE_END;
			tet.b_N[i] = SPACE_END;
			tet.c_N[i] = SPACE_END;
			tet.D[i] = SPACE_END;
		}
	}
	INTERSECT_TRI_P_INFO & operator =(const INTERSECT_TRI_P_INFO &src)
	{
		this ->tet = src.tet;
		this ->p = src.p;
		this ->t = src.t;
		return	*this;
	}
} INTERSECT_TRI_P_INFO;

typedef struct PointLight
{
	Point position;
	Point diffuseColor;
	float Kd; // == diffusePower
	Point specularColor;
	float Ks; // == specularPower
	float specularHardness;
	PointLight()
	{
		position = Point(SPACE_END, SPACE_END, SPACE_END);
		diffuseColor = Point(-1, -1, -1);
		Kd = 0;
		specularColor = Point(-1, -1, -1);
		Ks = 0;
		specularHardness = 0.0f;
	}
} PointLight;

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
	//vector<vector<vector<unsigned char> > > pixels2;
	//vector<vector<vector<unsigned char> > > pixels3;
	//vector<vector<vector<unsigned char> > > pixels4;

	//Sphere		*m_pSphere;
	list<SPHERE>	m_lSphere;
	int		m_nSphere;
	//Tetra
	list<TETRA>	m_lTetra;
	int		m_nTetra;
	list<TETRA_PP>	m_lTetraPP;
	int		m_nTetraPP;

	// Member Functions
	Raytracer(string nff_filename, string image_filename);
	~Raytracer();

	// Ray parser Reader Functions
	void Read();
	void Write();
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
	void DoPolyPP();

	HRESULT ImageAllocation(int nCol, int nRow);	//nCol == nWidth, nRow == nHeight
	COLOR3	RayIntersectWithObj(int i, int j);
	HRESULT	RaytracingInRange(int y_begin, int y_end);
	void	Raytracing();
	void RayIntersectWithSphere(double ray_origin[3], double ray_direction[3], SPHERE *sp,
			double *interP_x, double *interP_y, double *interP_z, double *_closest_t,
			double *normal_x, double *normal_y, double *normal_z, int ignoreID = -2, int nInside = -1);
	INTERSECT_TRI_INFO GetBarycentricCoordinate
			(double ray_origin[3], double ray_direction[3], int ignoreID = -2);
	INTERSECT_TRI_P_INFO GetBarycentricCoordinatePP
		(double ray_origin[3], double ray_direction[3], Point *normal, float *_alpha, float *_beta, float *_gamma, int ignoreID = -2);

	// Shadow Rays
	float *ShadowRays(Point &ray_origin, Point &normal, int nObj = -2, int ignoreID = -2);	
	COLOR3 GetPointLight(PointLight &light, Point &surfaceColor, Point &pos3D, Point &viewDir, Point &normal);
/// Check ignore ID part about TraceRAY......
	COLOR3 TraceRay(Point &rayOrigin, Point &rayDir, double _n1, int ray_depth, int _nHitObj = -1, int _ignoreID = -1, int nInside = -1);
	Point GetReflectRay(Point &normal, Point &rayDirection);
	Point GetRefractRay(Point &normal, Point &rayDirection, double _n1, double _n2);
};

Raytracer *g_ray;
static void * RaytracerThread(void *p)
{
	int *nThread = (int *)p;
	int upToY = (int)ceil(((float)(g_ray ->m_view.resolution.y) / g_nActiveThread));
	// Depending on the nThread number, go through different lines of y.
	g_ray ->RaytracingInRange(upToY * (*nThread), min(upToY * (*nThread+1), g_ray ->m_view.resolution.y));
	pthread_exit(NULL);
}

int main()
{
#ifdef PARAMETER_MODE
/*	if(argc != 4)
	{
		cout << "There must be 3 parameters" << endl <<
		"Usage>Raytracer.out <NFF file name> <PPM file name> <Thread No>" << endl;
		exit(-2);
	}

	// Object Oriented Programming style - Invoke Ray -> Load File.
	string arg1 = argv[1];
	string arg2 = argv[2];
	int no_threads = atoi(argv[3]);*/
string arg1 = "balls1.nff";//"balls1_transparent.nff";
string arg2 = "balls1.ppm";//"balls1_transparent.ppm";
int no_threads = 4;
#endif
	cout << "You indicated " << no_threads << " as # of threads." << endl;
	g_nActiveThread = no_threads;

	Raytracer ray(arg1, arg2);
	g_ray = &ray;
	ray.Read();
	ray.Raytracing();
	ray.Write();

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
	m_nTetraPP = 0;
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
		m_file >> sCommand;
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
				DoLightSource();
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
				if(sCommand == "pp")
					DoPolyPP();
				else
					DoPoly();
        		        break;
			case '\0' :
				break;
        		default:            /* unknown */
				cout << sCommand << endl;
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
	m_file >> m_view.camera.y;
	//cout << m_view.camera.y << endl;
	m_file >> m_view.camera.z;
	//cout << m_view.camera.z << endl;

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
	char s[256];
	m_file.getline(s, 255);
	char * pch;
	pch = strtok (s," ");
	float res[7];
	int i = 0;
	while (pch != NULL)
	{
		res[i] = atof(pch);
		pch = strtok(NULL, " ");
		i++;
	}
	m_view.lights[m_view.nLights].light.x = res[0];
	m_view.lights[m_view.nLights].light.y = res[1];
	m_view.lights[m_view.nLights].light.z = res[2];
	if(i > 4)	// light has colors
	{
		m_view.lights[m_view.nLights].lightColor.R = res[3];
		m_view.lights[m_view.nLights].lightColor.G = res[4];
		m_view.lights[m_view.nLights].lightColor.B = res[5];
		if(i == 7)
			m_view.lights[m_view.nLights].lightRadius = res[6];
		else
			m_view.lights[m_view.nLights].lightRadius = -1;
	}
	else
	{
		m_view.lights[m_view.nLights].lightColor.R = 1;
		m_view.lights[m_view.nLights].lightColor.G = 1;
		m_view.lights[m_view.nLights].lightColor.B = 1;
		m_view.lights[m_view.nLights].lightRadius = -1;
	}
	m_view.nLights++;
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
	float tmp;
	for(int i = 0 ; i < 8 ; i ++)
		m_file >> tmp;
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
	sp.ID = m_nSphere;
	sp.shd.color = m_shading.color;
	sp.shd.Kd = m_shading.Kd;
	sp.shd.Ks = m_shading.Ks;
	sp.shd.PhongPow = m_shading.PhongPow;
	sp.shd.Transmission = m_shading.Transmission;
	sp.shd.indexRefraction = m_shading.indexRefraction;

	m_lSphere.push_back(sp);
}

void Raytracer::DoPolyPP()
{
	int nVertices;
	m_file >> nVertices;
	if(nVertices < 3)
	{
		cout << "Polygons must have more than 3 vertices." << endl;
		exit(-6);
	}
	if(nVertices >= 3)
	{
		int nTriangles = nVertices - 3 + 1;
		float vertices[nVertices][3];
		float normal[nVertices][3];
		for(int i = 0 ; i < nVertices ; i++)
		{
			m_file >> vertices[i][0];
			m_file >> vertices[i][1];
			m_file >> vertices[i][2];
			m_file >> normal[i][0];
			m_file >> normal[i][1];
			m_file >> normal[i][2];
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
			float a[3], b[3], c[3], a_n[3], b_n[3], c_n[3];
			a[0] = vertices[TriangleVertexNo[i][0]][0];
			a_n[0] = normal[TriangleVertexNo[i][0]][0];
			a[1] = vertices[TriangleVertexNo[i][0]][1];
			a_n[1] = normal[TriangleVertexNo[i][0]][1];
			a[2] = vertices[TriangleVertexNo[i][0]][2];
			a_n[2] = normal[TriangleVertexNo[i][0]][2];
			b[0] = vertices[TriangleVertexNo[i][1]][0];
			b_n[0] = normal[TriangleVertexNo[i][1]][0];
			b[1] = vertices[TriangleVertexNo[i][1]][1];
			b_n[1] = normal[TriangleVertexNo[i][1]][1];
			b[2] = vertices[TriangleVertexNo[i][1]][2];
			b_n[2] = normal[TriangleVertexNo[i][1]][2];
			c[0] = vertices[TriangleVertexNo[i][2]][0];
			c_n[0] = normal[TriangleVertexNo[i][2]][0];
			c[1] = vertices[TriangleVertexNo[i][2]][1];
			c_n[1] = normal[TriangleVertexNo[i][2]][1];
			c[2] = vertices[TriangleVertexNo[i][2]][2];
			c_n[2] = normal[TriangleVertexNo[i][2]][2];
			// Now we have 3 vertices of the current triangle. and normal for each vertex, too!

			// Time to add the triangle.
			m_nTetraPP ++;
			TETRA_PP te_pp;
			te_pp.a[0] = a[0];
			te_pp.a[1] = a[1];
			te_pp.a[2] = a[2];
			te_pp.b[0] = b[0];
			te_pp.b[1] = b[1];
			te_pp.b[2] = b[2];
			te_pp.c[0] = c[0];
			te_pp.c[1] = c[1];
			te_pp.c[2] = c[2];
// //1.5 0 2.4 0.902861 0 0.429934
// if(b[0] > 1.4 && b[0] < 1.55 && b[2] < 2.41 && b[2] > 2.39
//  && b[1] < +0.01 && b[1] > -0.01)
// {
// cout <<"A" << endl; 
// cout << b_n[0] << " " << b_n[1] << " "<< b_n[2] << " " << endl;
// cout <<"B" << endl;
// }
			for(int tmpI = 0 ; tmpI < 3 ; tmpI ++)
			{
				te_pp.a_N[tmpI] = a_n[tmpI];
				te_pp.b_N[tmpI] = b_n[tmpI];
				te_pp.c_N[tmpI] = c_n[tmpI];
			}

		// Compute N and also D as in Ax + By + Cz + D = 0
		// I compute D because I think I can use the equation to process polygons with 4/5 vertices.
			Point polyNormal(       (te_pp.a_N[0] + te_pp.b_N[0] + te_pp.c_N[0]) / 3,
						(te_pp.a_N[1] + te_pp.b_N[1] + te_pp.c_N[1]) / 3,
						(te_pp.a_N[2] + te_pp.b_N[2] + te_pp.c_N[2]) / 3);
			CrossProduct(polyNormal.GetPoint(), te_pp.a, te_pp.D);
			te_pp.D[0] = -te_pp.D[0];
			te_pp.D[1] = -te_pp.D[1];
			te_pp.D[2] = -te_pp.D[2];
			te_pp.ID = m_nTetraPP;
			te_pp.shd = m_shading;
			m_lTetraPP.push_back(te_pp);
		}
	}
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
			double b_a[3], c_a[3];
			MatrixSubtract1X3(te.b, te.a, b_a);
			MatrixSubtract1X3(te.c, te.a, c_a);
			CrossProduct(b_a, c_a, te.N);
			CrossProduct(te.N, te.a, te.D);
			te.D[0] = -te.D[0];
			te.D[1] = -te.D[1];
			te.D[2] = -te.D[2];
			te.ID = m_nTetra;
			te.shd = m_shading;
			m_lTetra.push_back(te);
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

// return distance .. if ray is blocked, return -1;
float *Raytracer::ShadowRays(Point &ray_origin, Point &normal, int nObj, int ignoreID)
{
	float *array_t = new float[m_view.nLights];
	for(int e = 0 ; e < m_view.nLights ; e++)
		array_t[e] = -1;
	for(int i = 0 ; i < m_view.nLights ; i++)
	{
	// If the ray goes backwards to the other side, you should not even shoot a ray.
	// because light cannot penetrate to reach light behind the object.
		Point ray(m_view.lights[i].light);
		Point ray_direction(ray - ray_origin);
		double distanceToLight = ray_direction.GetLengthVector();
		Point normalizedDir[9];
		normalizedDir[0] = ray_direction;
		normalizedDir[0].NormalizeVector();
		float bBackCull = Point::DotProduct(normalizedDir[0], normal);
		if(bBackCull < 0.0f)
		{
			array_t[i] = -1;
			continue;
		}

		Point *tmpup = new Point(0,1,0);
		Point *invec = new Point;
		Point::CrossProduct(*tmpup, normalizedDir[0], *invec);
		invec ->NormalizeVector();
		Point *invec2 = new Point;
		Point::CrossProduct(*invec, normalizedDir[0], *invec2);
		invec2 ->NormalizeVector();
		Point *invec3 = new Point;
		Point *invec4 = new Point;
Point *invec5 = new Point;
Point *invec6 = new Point;
Point *invec7 = new Point;
Point *invec8 = new Point;
		double lightR;
		if(m_view.lights[i].lightRadius < 0)
			lightR = 0.01;
		else
			lightR = m_view.lights[i].lightRadius;
		*invec = *invec * lightR;
		*invec2 = *invec2 * lightR;
		*invec3 = *invec * -1;
		*invec4 = *invec2 * -1;
*invec5 = *invec * (lightR * 0.5);
*invec6 = *invec2 * (lightR * 0.5);
*invec7 = *invec5 * -1;
*invec8 = *invec6 * -1;

		*invec = *invec + m_view.lights[i].light;
		*invec2 = *invec2 + m_view.lights[i].light;
		*invec3 = *invec3 + m_view.lights[i].light;
		*invec4 = *invec4 + m_view.lights[i].light;
*invec5 = *invec5 + m_view.lights[i].light;
*invec6 = *invec6 + m_view.lights[i].light;
*invec7 = *invec7 + m_view.lights[i].light;
*invec8 = *invec8 + m_view.lights[i].light;

		normalizedDir[1] = *invec - ray_origin; normalizedDir[2] = *invec2 - ray_origin;
		normalizedDir[3] = *invec3 - ray_origin; normalizedDir[4] = *invec4 - ray_origin;
		normalizedDir[5] = *invec5 - ray_origin; normalizedDir[6] = *invec6 - ray_origin;
		normalizedDir[7] = *invec7 - ray_origin; normalizedDir[8] = *invec8 - ray_origin;
		delete invec, invec2, invec3, invec4, invec5, invec6, invec7, invec8;
		INTERSECT_SPHERE_INFO sp_info[9];
		INTERSECT_TRI_INFO tri_info[9];
		INTERSECT_TRI_P_INFO tri_p_info[9];

		int sp_id, tri_id, tri_p_id;
		if(nObj == 1)
			sp_id = ignoreID;
		else if(nObj == 2)
			tri_id = ignoreID;
		else if(nObj == 3)
			tri_p_id = ignoreID;
		float *alpha = new float;
		float *beta = new float;
		float *gamma = new float;
		Point *normalBarycentric = new Point;
		double normal_x, normal_y, normal_z;
		for(int e = 0 ; e < 9 ; e ++)
		{
			RayIntersectWithSphere(ray_origin.GetPoint(), normalizedDir[e].GetPoint(), &sp_info[e].sp,
				&sp_info[e].p.x, &sp_info[e].p.y, &sp_info[e].p.z, &sp_info[e].t, &normal_x, &normal_y, &normal_z, sp_id);
			tri_info[e] = GetBarycentricCoordinate(ray_origin.GetPoint(), normalizedDir[e].GetPoint(), tri_id);
			tri_p_info[e] = GetBarycentricCoordinatePP
				(ray_origin.GetPoint(), normalizedDir[e].GetPoint(), normalBarycentric, alpha, beta, gamma, tri_p_id);
		}
		delete alpha; delete beta; delete gamma; delete normalBarycentric;
		array_t[i] = 0;
		for(int e = 0 ; e < 9 ; e++)
		{
			if(sp_info[e].t < distanceToLight || tri_info[e].t < distanceToLight || tri_p_info[e].t < distanceToLight) // Blocked
				;
			else
				array_t[i] += 1.0/9;
		}
	}
	return array_t;
}

void Raytracer::RayIntersectWithSphere(double ray_origin[3], double ray_direction[3], SPHERE *sp, double *interP_x, double *interP_y,
				double *interP_z, double *_closest_t, double *normal_x, double *normal_y,
				double *normal_z, int ignoreID, int nInside)
{
	//  Now we have Us, Vs, S, so just compute t right away.
	// Compute parameter t for all the spheres in LIST.
	// and then get the closest t. We are going to use the closest t to display colors.
	double s[3];
	Point pS;
	pS = ray_direction;
	pS += ray_origin;
	s[0] = pS.x; s[1] = pS.y; s[2] = pS.z;
	double d[3]; // == ray direction
	double e[3]; // == eye(ray origin)
	e[0] = ray_origin[0] ; e[1] = ray_origin[1]; e[2] = ray_origin[2];
	// d = s - e;
	d[0] = s[0] - e[0]; d[1] = s[1] - e[1]; d[2] = s[2] - e[2];
	double d_dot_d;
	DotProduct(d, d, &d_dot_d);
	double minus_d[3];
	minus_d[0] = -d[0]; minus_d[1] = -d[1]; minus_d[2] = -d[2];
	// e_c varies for each sphere. so I must compute e_c for all the spheres.
	// Radius also varies for each sphere.
	// Conclusion ! > I must compute R, e_c for each sphere.
	double e_c[3];
	double radius;
	double c[3];
	double closest_t = SPACE_END;	// closest t
//	float closest_c[3] = {SPACE_END, SPACE_END, SPACE_END}; // center
//	float closest_R = SPACE_END;	// radius of the closest sphere.
	Point computedNormal;
	INTERSECT_SPHERE_INFO sp_info;	// Automatic initialization.
	for(list<SPHERE>::iterator it = m_lSphere.begin() ; it != m_lSphere.end() ; it++)
	{
	/// In case of refraction / reflection ray, you must change the code. ==> Change it so that it will not ignore the backface sphere.
	/// because it should do refraction when it gets out of the sphere in which the ray is.
	/// For example if((*it).ID == ignoreID && Dot(interP.surfaceNormal, rayDirection) > 0)
		if(nInside == -1) // Normal Case
		{
			if((*it).ID == ignoreID)
			{
				continue;
			}
		}
		else if(nInside == 0) // Refracted Ray But it's outside.
		{
		}
		else	// Refracted Ray but it's inside.
		{	// Hit the backside of the sphere, but when it hits, it will not be reflected but refracted again.
		}
		// c = center of each sphere
		c[0] = (*it).x;	c[1] = (*it).y; c[2] = (*it).z;
		radius = (*it).radius;
		e_c[0] = e[0] - c[0]; e_c[1] = e[1] - c[1]; e_c[2] = e[2] - c[2];
		double determinant = 0.0f;
		double tmp1;
		DotProduct(d, e_c, &tmp1);
		tmp1 = tmp1 * tmp1;
		double tmp2;
		DotProduct(e_c, e_c, &tmp2);
		tmp2 -= (radius * radius);
		double tmp3;
		tmp3 = d_dot_d * tmp2;
		determinant = tmp1 - tmp3;
		// Just for debugging. I am wary of floating point precision issues because it might not be 0 but 0.0000001
//		if(0.0000001 > determinant && determinant > -0.0000001)
//		{
//			cout << "There might be precision issue. Check line 463" << endl;
//		}
		if(determinant > 0)	// 2 roots
		{
///				if((*it).ID == ignoreID && Point::DotProduct(ray_direction,computedNormal) > 0)
///					continue;
			double t1, t2;
			double term1;
			DotProduct(minus_d, e_c, &term1);
			t1 = term1 + sqrt(determinant);
			t1 /= d_dot_d;
			t2 = term1 - sqrt(determinant);
			t2 /= d_dot_d;
			if(t1 < 0)
				t1 = SPACE_END;
			if(t2 < 0)
				t2 = SPACE_END;
			if(t2 > t1 && closest_t > t1)
			{	// Whenever we find the t, compute N and save it in the global variable.
				Point center = c;
				computedNormal = ray_direction;
				computedNormal *= t1;
				computedNormal += ray_origin;
				computedNormal -= center;
				computedNormal.NormalizeVector();
//				if((*it).ID == ignoreID && Point::DotProduct(ray_direction,computedNormal) > 0)
//					continue;
				sp_info.t = t1;
				sp_info.sp = (*it);
				sp_info.sp.x = c[0];
				sp_info.sp.y = c[1];
				sp_info.sp.z = c[2];
				sp_info.sp.radius = radius;
				sp_info.sp.ID = (*it).ID;
				sp_info.t = closest_t = t1;
			}
			//closest_t = SMALLER(closest_t, t1);
			if(t1 > t2 && closest_t > t2)
			{
				Point center = c;
				computedNormal = ray_direction;
				computedNormal *= t2;
				computedNormal += ray_origin;
				computedNormal -= center;
				computedNormal.NormalizeVector();
				sp_info.t = t2;
				sp_info.sp = (*it);
				sp_info.sp.x = c[0];
				sp_info.sp.y = c[1];
				sp_info.sp.z = c[2];
				sp_info.sp.radius = radius;
				sp_info.sp.ID = (*it).ID;
				sp_info.t = closest_t = t2;
			}
			//closest_t = SMALLER(closest_t, t2);
		}
		else if(determinant == 0)	// only 1 root
		{
			double t;
			DotProduct(minus_d, e_c, &t);
			t /= d_dot_d;
			if(closest_t > t && t > 0)
			{
				Point center = c;
				computedNormal = ray_direction;
				computedNormal *= t;
				computedNormal += ray_origin;
				computedNormal -= center;
				computedNormal.NormalizeVector();
				sp_info.t = t;
				sp_info.sp = (*it);
				sp_info.sp.x = c[0];
				sp_info.sp.y = c[1];
				sp_info.sp.z = c[2];
				sp_info.sp.radius = radius;
				sp_info.sp.ID = (*it).ID;
				sp_info.t = closest_t = t;
			}
			//closest_t = SMALLER(closest_t, t);
		}
		else				// no roots -- Not intersects!
		{
			;
		}
	}
	*normal_x = computedNormal.x; *normal_y = computedNormal.y; *normal_z = computedNormal.z;
	Point pDir(ray_direction), pOri(ray_origin);
	sp_info.p = pOri + pDir * closest_t;
	sp ->x = sp_info.sp.x; sp ->y = sp_info.sp.y; sp ->z = sp_info.sp.z;
	sp ->radius = sp_info.sp.radius; sp ->ID = sp_info.sp.ID; sp ->shd = sp_info.sp.shd;
	*interP_x = sp_info.p.x; *interP_y = sp_info.p.y; *interP_z = sp_info.p.z;
	*_closest_t = closest_t;

	return;
}

COLOR3 Raytracer::RayIntersectWithObj(int i, int j)
{
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
	double w[3];
	double fLookAt[3], fCameraPos[3];
	fLookAt[0] = m_view.lookAt.x;
	fLookAt[1] = m_view.lookAt.y;
	fLookAt[2] = m_view.lookAt.z;
	fCameraPos[0] = m_view.camera.x;
	fCameraPos[1] = m_view.camera.y;
	fCameraPos[2] = m_view.camera.z;
	MatrixSubtract1X3(fLookAt, fCameraPos, w);//w = gaze(viewing) vector
	NormalizeVector(w, w);
	// Camera up vector.. get it! - camera_up[3]
	// To get it, I set an arbitrary vector and then get the other 2 vectors from it.
	double arbitrary_vector[3]= {m_view.cameraUp.x, m_view.cameraUp.y, m_view.cameraUp.z};
	double camera_right[3];
	CrossProduct(w, arbitrary_vector, camera_right);
	NormalizeVector(camera_right, camera_right);
	double camera_up[3];
	CrossProduct(camera_right, w, camera_up);
	NormalizeVector(camera_up, camera_up);

	// Now we know camera_up, w(viewing direction), camera_right)
	// s = e + Us*U + Vs*V - d*W (s = origin(e) + direction(Us*U + Vs*V - d*w))
	double ray_direction[3] = {0.0f, 0.0f, 0.0f};// = {camera_right, camera_up, w};
	camera_right[0] *= Us; camera_right[1] *= Us; camera_right[2] *= Us;
	// Image plane starts from top, which is (0,0) to bottom(512,512)
	camera_up[0] *= -Vs; camera_up[1] *= -Vs; camera_up[2] *= -Vs;
	MatrixSum1X3(camera_right, camera_up, ray_direction);
	w[0] *= m_view.dNear; w[1] *= m_view.dNear; w[2] *= m_view.dNear;
	MatrixSum1X3(ray_direction, w, ray_direction);
	NormalizeVector(camera_right, camera_right);
	NormalizeVector(camera_up, camera_up);
	NormalizeVector(w, w);
	NormalizeVector(ray_direction, ray_direction);
	double ray_origin[3];
	ray_origin[0] = m_view.camera.x; ray_origin[1] = m_view.camera.y; ray_origin[2] = m_view.camera.z;

	Point pRayOrigin(ray_origin);
	Point pRayDir(ray_direction);
	double closest_t = SPACE_END;
	INTERSECT_SPHERE_INFO sp_info;
	INTERSECT_TRI_INFO tri_info;
	INTERSECT_TRI_P_INFO tri_p_info;
	float *alpha = new float; *alpha = 0;
	float *beta = new float; *beta = 0;
	float *gamma = new float; *gamma = 0;
	Point *normal_from_Barycentric = new Point;
	double normal_x, normal_y, normal_z;
	RayIntersectWithSphere(ray_origin, ray_direction, &sp_info.sp, &sp_info.p.x,
		&sp_info.p.y, &sp_info.p.z, &sp_info.t, &normal_x, &normal_y, &normal_z);
	tri_info = GetBarycentricCoordinate(ray_origin, ray_direction);
	tri_p_info = GetBarycentricCoordinatePP(ray_origin, ray_direction, normal_from_Barycentric, alpha, beta, gamma);

// HERE If intersected, shoot a ray onto the light sources. and create shadow first.
	//shoot a ray onto light sources.

	int nHitObj = 0; // nHitObj == 1(Sphere) , nHitObj == 2(Triangle) , nHitObj == 3 (Triangle with PP)
	Point interP;
	if(sp_info.t < tri_info.t && sp_info.t < tri_p_info.t)
	{
		nHitObj = 1;
		closest_t = sp_info.t;
		interP = sp_info.p;
	}
	else if(tri_info.t < sp_info.t && tri_info.t < tri_p_info.t)
	{
		nHitObj = 2;
		closest_t = tri_info.t;
		interP = tri_info.p;
	}
	else if(tri_p_info.t < sp_info.t && tri_p_info.t < tri_info.t)
	{
		nHitObj = 3;
		closest_t = tri_p_info.t;
		interP = tri_p_info.p;
	}

	// At this point, we know the closest t. If t == 100000.0(MACRO:SPACE_END), the ray does NOT intersect with any spheres.
	if(closest_t == SPACE_END)
	{
		// Never intersect with any spheres that we have.
		// In this case, just return color = {-1, -1, -1};
		COLOR3 color(-1, -1, -1);
		return color;
	}
	else if(closest_t != SPACE_END)
	{
		// I am going to generate a ray and shoot it at each light source.
		COLOR3 color;
//		float Kd, Ks, PhongPow;
		Point interP;
		Point normal;
		SHADING shader;
		int ignoreID = -2;
		if(nHitObj == 1)
		{
			interP = sp_info.p;
			color = sp_info.sp.shd.color;
			Point center(sp_info.sp.x, sp_info.sp.y, sp_info.sp.z);
			Point tmp_normal;
			tmp_normal = interP - center;
			tmp_normal.NormalizeVector();
			//normal = tmp_normal;
			normal = Point(normal_x, normal_y, normal_z);
			ignoreID = sp_info.sp.ID;
			shader = sp_info.sp.shd;
		}
		else if(nHitObj == 2)
		{
			interP = tri_info.p;
			color = tri_info.tet.shd.color;
			normal = tri_info.tet.N;
			ignoreID = tri_info.tet.ID;
			shader = tri_info.tet.shd;
		}
		else if(nHitObj == 3)
		{
//cout << *gamma << " " << *beta << " " << 1 - *beta - *gamma << endl;
			interP = tri_p_info.p;
			color = tri_p_info.tet.shd.color;
			// kind of complicated to get the normal of surface(polygon.). Get the normals of each vertex and get avg.
			normal = *normal_from_Barycentric;
			normal.NormalizeVector();
			ignoreID = tri_p_info.tet.ID;
			shader = tri_p_info.tet.shd;
		}
		Point ray_origin(interP);
		//ray direction is light sources.
		float *pT = ShadowRays(interP, normal, nHitObj, ignoreID);
		COLOR3 wholesumColor(0,0,0);
		Point surfaceColor(color);
		for(int i = 0 ; i < m_view.nLights ; i++)
		{
			if(pT[i] > 0)
			{
				PointLight light0;
				light0.position = m_view.lights[i].light;
				light0.diffuseColor = m_view.lights[i].lightColor;
				light0.Kd = shader.Kd;	// get it from the object surface.
				light0.specularColor = m_view.lights[i].lightColor;
				light0.Ks = shader.Ks;
				light0.specularHardness = shader.PhongPow;
				// We have to consider the surface color on the intersecting point.

				Point pCamera(m_view.camera);
				Point viewDir = pCamera - interP;
				COLOR3 pointColor(0,0,0);
				pointColor = GetPointLight(light0, surfaceColor, interP, viewDir, normal);
				wholesumColor.R += pointColor.R * pT[i];
				wholesumColor.G += pointColor.G * pT[i];
				wholesumColor.B += pointColor.B * pT[i];
			}
		}
		double division = sqrt(m_view.nLights);
		division = 1.0 / division;
		wholesumColor.R *= division; wholesumColor.G *= division; wholesumColor.B *= division;

		COLOR3 reflect_color(0,0,0);
		Point normalizedDir = pRayDir;
		normalizedDir.NormalizeVector();
		normal.NormalizeVector();
		if(shader.Ks > 0 && Point::DotProduct(normalizedDir, normal) < 0)
		{
			double air_n = 1.0;// current refract index; I assume that the air refract index is 1.0.
			Point reflectedRay = GetReflectRay(normal, pRayDir);
			reflect_color = TraceRay(interP, reflectedRay, air_n, 1, nHitObj, ignoreID);	// Reflection --> nInside = -1
		}
		COLOR3 refract_color(0,0,0);
		if(shader.Transmission > 0)
		{
			double air_n = 1.0;
			Point refractedRay = GetRefractRay(normal, normalizedDir, air_n, shader.indexRefraction);
			if(refractedRay.x != SPACE_END || refractedRay.y != SPACE_END || refractedRay.z != SPACE_END)
				// No refraction if refractedRay sqrt returns nan.
				refract_color = TraceRay(interP, refractedRay, air_n, 5, nHitObj, ignoreID, 1); // nInside == 1
		}
		wholesumColor.R = wholesumColor.R + reflect_color.R * shader.Ks + refract_color.R * shader.Transmission/2;
		wholesumColor.G = wholesumColor.G + reflect_color.G * shader.Ks + refract_color.G * shader.Transmission/2;
		wholesumColor.B = wholesumColor.B + reflect_color.B * shader.Ks + refract_color.B * shader.Transmission/2;
/*Point a = reflect_color;
Point b = refract_color;
if(a.IsZero() == 0)
     a.Print();
if(b.IsZero() == 0)
b.Print();*/
		delete pT;
		delete alpha; delete beta;delete gamma; delete normal_from_Barycentric;
		return	wholesumColor;
	}

	COLOR3 color(-1 , -1 , -1);
	delete alpha; delete beta;delete gamma; delete normal_from_Barycentric;
	return	color;
}


void Raytracer::Raytracing()
{
	pthread_t thread[g_nActiveThread];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	int rc;

	g_nActiveThread = min(g_nActiveThread, m_view.resolution.y);
	ImageAllocation(m_view.resolution.x, m_view.resolution.y);
	int nThreadID[g_nActiveThread];
	for(int i = 0 ; i < g_nActiveThread ; i++)
		nThreadID[i] = i;
	struct timeval tim1, tim2;
	gettimeofday(&tim1, NULL);
	for(int nThread = 0 ; nThread < g_nActiveThread ; nThread++)
	{
		rc = pthread_create(&thread[nThread], &attr, RaytracerThread, (void*) &nThreadID[nThread]);
		if(rc != 0)
		{
			cout << "Thread Creation Failed. Error Code : " << rc << endl;
			exit(-9);
		}
	}
	int *nStatus;
	for(int nThread = 0 ; nThread < g_nActiveThread ; nThread ++)
	{
		rc = pthread_join(thread[nThread], (void **)&nStatus);
		if(rc != 0)
		{
			cout << "Thread Join Failed. Error Code : " << rc << endl;
			exit(-8);
		}
	}
	gettimeofday(&tim2, NULL);
	long mtime, seconds , useconds;
	seconds = tim2.tv_sec - tim1.tv_sec;
	useconds = tim2.tv_usec - tim1.tv_usec;
	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
	cout << "Elapsed time : " << mtime << endl;
}


HRESULT	Raytracer::RaytracingInRange(int y_begin, int y_end)
{
//	Here, Use thread or whatever.. Just traverse with i, j.
//	i = [l, r], r = [t,b]
	// Now save the data into the vectors(== array)
	int resX = m_view.resolution.x;
	for(int j = y_begin ; j < y_end ; j++)
	{
		for(int i = 0 ; i < resX ; i++)
		{
			COLOR3 color(0,0,0);
			color = RayIntersectWithObj(i, j);
			if(color.R == -1 && color.G == -1 && color.B == -1)
			{
				// Never intersect
				// Later, in this part, color should be the background color.
				color.R = m_background.R;
				color.G = m_background.G;
				color.B = m_background.B;
			}
			color.R *= 255; color.G *= 255; color.B *= 255;
			CLIPPING(color.R);
			CLIPPING(color.G);
			CLIPPING(color.B);
			unsigned char R = (unsigned char)color.R;
			unsigned char G = (unsigned char)color.G;
			unsigned char B = (unsigned char)color.B;
			pixels[j][i][0] = R;
			pixels[j][i][1] = G;
			pixels[j][i][2] = B;
		}
	}
}


void MatrixMult3X3(float M1[][3], float M2[][3], float resultM[][3])
{
	int i = 0, j = 0 , k = 0;
	for(i = 0 ; i < 3 ; i++)
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

void DotProduct(double a[3], double b[3], double *result)
{
	*result = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void MatrixSubtract1X3(float M1[3], float M2[3], float resultM[3])
{
	resultM[0] = M1[0] - M2[0];
	resultM[1] = M1[1] - M2[1];
	resultM[2] = M1[2] - M2[2];
}

void MatrixSubtract1X3(double M1[3], double M2[3], double resultM[3])
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

void MatrixSum1X3(double M1[3], double M2[3], double resultM[3])
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
void NormalizeVector(double M[3], double resultM[3])
{
	double length = M[0] * M[0] + M[1] * M[1] + M[2] * M[2];
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

void CrossProduct(double M1[3], double M2[3], double resultM[3])
{
	resultM[0] = M1[1] * M2[2] - M1[2] * M2[1];
	resultM[1] = M1[2] * M2[0] - M1[0] * M2[2];
	resultM[2] = M1[0] * M2[1] - M1[1] * M2[0];
}

INTERSECT_TRI_INFO Raytracer::GetBarycentricCoordinate
			(double ray_origin[3], double ray_direction[3], int ignoreID)
{
	INTERSECT_TRI_INFO tri_info;
	double closest_t = SPACE_END;
	for(list<TETRA>::iterator it = m_lTetra.begin() ; it != m_lTetra.end() ; it++)
	{
		TETRA tet = (*it);
		if(tet.ID == ignoreID)
			continue;
		// Compute A
		float a, b, c, d, e, f, g, h, i;
		a = tet.a[0] - tet.b[0];
		b = tet.a[1] - tet.b[1];
		c = tet.a[2] - tet.b[2];
		d = tet.a[0] - tet.c[0];
		e = tet.a[1] - tet.c[1];
		f = tet.a[2] - tet.c[2];
		g = ray_direction[0];
		h = ray_direction[1];
		i = ray_direction[2];
		float j, k, l;
		j = tet.a[0] - ray_origin[0];
		k = tet.a[1] - ray_origin[1];
		l = tet.a[2] - ray_origin[2];

		// Compute M
		double M = a * (e * i - h *f) + b * (g * f - d * i) + c * (d * h - e * g);
		// if M == 0 ==> no solution ==> not intersect
		if(M == 0)
		{
			cout << "The triangle is degenerate or parallel to the ray." << endl;
			continue;
			//return	-1;
		}
		// if M != 0 ==> get alpha, beta, t
		double t = f * (a * k - j * b) + e * (j * c - a * l) + d * (b * l - k * c);
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

		if(closest_t > t && t > 0)
		{
			tri_info.tet = tet;
			tri_info.t = closest_t = t;
			Point pRayOrigin(ray_origin);
			Point pRayDir(ray_direction);
			Point pa = tet.a;
			Point pb = tet.b;
			Point pc = tet.c;
			Point tmp = pa * (1-gamma-beta);
			Point tmp2 = pb * beta;
			Point tmp3 = pc * gamma;
			tri_info.p = tmp;//pRayOrigin + pRayDir * closest_t;
			tri_info.p += tmp2;
			tri_info.p += tmp3;
		}
	}
// calculate intersecting point p
	return	tri_info;
}

INTERSECT_TRI_P_INFO Raytracer::GetBarycentricCoordinatePP
		(double ray_origin[3], double ray_direction[3], Point *normal, float *_alpha, float *_beta, float *_gamma, int ignoreID)
{
	double closest_t = SPACE_END;
	Point p;
	TETRA_PP storageTet;
	for(list<TETRA_PP>::iterator it = m_lTetraPP.begin() ; it != m_lTetraPP.end() ; it++)
	{
		if((*it).ID == ignoreID)
			continue;
		// Compute A
		TETRA_PP tet = (*it);
		double a, b, c, d, e, f, g, h, i;
		a = tet.a[0] - tet.b[0];
		b = tet.a[1] - tet.b[1];
		c = tet.a[2] - tet.b[2];
		d = tet.a[0] - tet.c[0];
		e = tet.a[1] - tet.c[1];
		f = tet.a[2] - tet.c[2];
		g = ray_direction[0];
		h = ray_direction[1];
		i = ray_direction[2];
		double j, k, l;
		j = tet.a[0] - ray_origin[0];
		k = tet.a[1] - ray_origin[1];
		l = tet.a[2] - ray_origin[2];

		// Compute M
		double M = a * (e * i - h *f) + b * (g * f - d * i) + c * (d * h - e * g);
		// if M == 0 ==> no solution ==> not intersect
		if(M == 0)
		{
			cout << "The triangle is degenerate or parallel to the ray." << endl;
			continue;
			//return	-1;
		}
		// if M != 0 ==> get alpha, beta, t
		double t = f * (a * k - j * b) + e * (j * c - a * l) + d * (b * l - k * c);
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
		*_alpha = 1 - beta - gamma;
		*_beta = beta;
		*_gamma = gamma;

		if(closest_t > t && t > 0)
		{
			storageTet = tet;
			closest_t = t;
			Point pa = tet.a;
			Point pb = tet.b;
			Point pc = tet.c;
			Point tmp1 = pa * (1-beta-gamma);
			Point tmp2 = pb * beta;
			Point tmp3 = pc * gamma;
			p = tmp1 + tmp2 + tmp3;//pRayOrigin + pRayDir * closest_t;
			Point a_n = tet.a_N; Point b_n = tet.b_N; Point c_n = tet.c_N;
			*normal = a_n * (1-beta-gamma) + b_n* beta + c_n * gamma;
		}
	}

	INTERSECT_TRI_P_INFO tri_info;
	tri_info.tet = storageTet;
	tri_info.t = closest_t;
	tri_info.p = p;

	return	tri_info;
}

void Raytracer::Write()
{
	m_outFile.open(m_output_image_filename.c_str());
	int resY, resX;
	resX = m_view.resolution.x;
	resY = m_view.resolution.y;
	unsigned char buffer[resX * 3];
	m_outFile << "P6\n" << resX << " " << resY << "\n255\n";
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

COLOR3 Raytracer::GetPointLight(PointLight &light, Point &surfaceColor, Point &pos3D, Point &viewDir, Point &normal)
{
	Point DiffuseColor;
	Point SpecularColor;
	normal.NormalizeVector();
	Point normalizedViewDir = viewDir;
	normalizedViewDir *= -1;
	normalizedViewDir.NormalizeVector();
	if(light.Kd > 0)
	{
		Point lightDir = light.position - pos3D; // FIND THE VECTOR BETWEEN THE 3D POSITION IN SPACE OF THE SURFACE
		float distance = lightDir.GetLengthVector(); // GET THE DISTANCE OF THIS VECTOR
//		distance = distance * distance; //+ 2 * distance + 0.5; // USES INVERSE SQUARE FOR DISTANCE ATTENUATION
		lightDir.NormalizeVector(); // NORMALIZE THE VECTOR

		// INTENSITY OF THE DIFFUSE LIGHT	
                // KEEP WITHIN THE 0-1 RANGE
                // DOT PRODUCT OF THE LIGHT DIRECTION VECTOR AND THE SURFACE NORMAL
		float i = Point::DotProduct(lightDir, normal);
		i=(i>1.0)?1.0:i;
		i=(i<0.0)?0.0:i;

		// CALCULATE THE DIFFUSE LIGHT FACTORING IN LIGHT COLOUR, POWER AND THE ATTENUATION
		double distance_factor = 0.0034868992* distance * distance - 0.03794031 * distance + 0.9531702318;
// 0.03 --> For balls // * 0.0035; --> For Teapot
		
		DiffuseColor = light.diffuseColor * surfaceColor * (i * light.Kd / distance_factor); //Later divide by distance.

                //CALCULATE THE HALF VECTOR BETWEEN THE LIGHT VECTOR AND THE VIEW VECTOR.
		//THIS IS CHEAPER THEN CALCULATING THE ACTUAL REFLECTIVE VECTOR
                Point h = lightDir + normalizedViewDir;//normalize(lightDir + viewDir);
		h.NormalizeVector();

	        // INTENSITY OF THE SPECULAR LIGHT	
                // DOT PRODUCT OF NORMAL VECTOR AND THE HALF VECTOR TO THE POWER OF THE SPECULAR HARDNESS
		float specularTerm = Point::DotProduct(normal, h);
		specularTerm=(specularTerm>1.0)?1.0:specularTerm;
		specularTerm=(specularTerm<0.0)?0.0:specularTerm;
		i = pow(specularTerm, light.specularHardness);
// Point minus_view = viewDir * -1;	
// if(Point::SimilarityCheck(normal, minus_view) == 1)
// cout << "A" << endl;

		// CALCULATE THE SPECULAR LIGHT FACTORING IN LIGHT SPECULAR COLOUR, POWER AND THE ATTENUATION
		SpecularColor = light.diffuseColor * (i * light.Ks / distance_factor); // / distance;
	}
	Point color = DiffuseColor + SpecularColor;
	COLOR3 out;
	out.R = color.x; out.G = color.y; out.B = color.z;
	return out;
}

COLOR3 Raytracer::TraceRay(Point &rayOrigin, Point &rayDir, double _n1, int ray_depth, int _nHitObj, int _ignoreID, int nInside)
{
	if(ray_depth == 0)
	{
		COLOR3 clr(0,0,0);
		return	clr;
	}

	INTERSECT_SPHERE_INFO sp_info;
	INTERSECT_TRI_INFO tri_info;
	INTERSECT_TRI_P_INFO tri_p_info;
	float *alpha = new float; *alpha = 0;
	float *beta = new float; *beta = 0;
	float *gamma = new float; *gamma = 0;
	Point *normal_from_Barycentric = new Point;
	double normal_x, normal_y, normal_z;
	double ray_ori[3] = {rayOrigin.x, rayOrigin.y, rayOrigin.z};
	double ray_dir[3] = {rayDir.x, rayDir.y, rayDir.z};
	if(_nHitObj == 1)
		RayIntersectWithSphere(ray_ori, ray_dir, &sp_info.sp,
			&sp_info.p.x, &sp_info.p.y, &sp_info.p.z, &sp_info.t,
			&normal_x, &normal_y, &normal_z, _ignoreID, nInside);
	else
		RayIntersectWithSphere(ray_ori, ray_dir, &sp_info.sp,
			&sp_info.p.x, &sp_info.p.y, &sp_info.p.z, &sp_info.t, &normal_x, &normal_y, &normal_z);
	if(_nHitObj == 2)
		tri_info = GetBarycentricCoordinate(ray_ori, ray_dir, _ignoreID);
	else
		tri_info = GetBarycentricCoordinate(ray_ori, ray_dir);
	if(_nHitObj == 3)
		tri_p_info = GetBarycentricCoordinatePP
			(ray_ori, ray_dir, normal_from_Barycentric, alpha, beta, gamma, _ignoreID);
	else
		tri_p_info = GetBarycentricCoordinatePP
			(ray_ori, ray_dir, normal_from_Barycentric, alpha, beta, gamma);

	int nHitObj = 0; // nHitObj == 1(Sphere) , nHitObj == 2(Triangle) , nHitObj == 3 (Triangle with PP)
	int nIgnoreID = -1;
	Point interP;
	Point normal;
	double closest_t;
	SHADING shader;
	if(sp_info.t == SPACE_END && tri_info.t == SPACE_END && tri_p_info.t == SPACE_END)
	{
		COLOR3 clr_background;
		clr_background = m_background;
		return COLOR3(0,0,0);
	}
	if(sp_info.t < tri_info.t && sp_info.t < tri_p_info.t)
	{
		nHitObj = 1;
		nIgnoreID = sp_info.sp.ID;
		closest_t = sp_info.t;
		interP = sp_info.p;
		shader = sp_info.sp.shd;
		normal = Point(normal_x, normal_y, normal_z);
if(nInside == 1)
	nInside = 2;	// Special case - It's about to get out of the sphere from the inside of it. --> In this case, it does not get color but only refract.(No reflection either.)
	}
	else if(tri_info.t < sp_info.t && tri_info.t < tri_p_info.t)
	{
		nHitObj = 2;
		nIgnoreID = tri_info.tet.ID;
		closest_t = tri_info.t;
		interP = tri_info.p;
		shader = tri_info.tet.shd;
		normal = tri_info.tet.N;
	}
	else if(tri_p_info.t < sp_info.t && tri_p_info.t < tri_info.t)
	{
		nHitObj = 3;
		nIgnoreID = tri_p_info.tet.ID;
		closest_t = tri_p_info.t;
		interP = tri_p_info.p;
		shader = tri_p_info.tet.shd;
		normal = *normal_from_Barycentric;
	}

	Point shadow_ray_origin(interP);
	//ray direction is light sources.
	float *pT = ShadowRays(shadow_ray_origin, normal, nHitObj, nIgnoreID);
	COLOR3 wholesumColor(0,0,0);
	Point surfaceColor(shader.color);
	for(int i = 0 ; i < m_view.nLights && nInside != 2 ; i++)
	{
		if(pT[i] > 0)
		{
			PointLight light0;
			light0.position = m_view.lights[i].light;
			light0.diffuseColor = m_view.lights[i].lightColor;
			light0.Kd = shader.Kd;	// get it from the object surface.
			light0.specularColor = m_view.lights[i].lightColor;
			light0.Ks = shader.Ks;
			light0.specularHardness = shader.PhongPow;
			// We have to consider the surface color on the intersecting point.

			Point pCamera(m_view.camera);
			Point viewDir = pCamera - interP;
			COLOR3 pointColor(0,0,0);
			pointColor = GetPointLight(light0, surfaceColor, interP, viewDir, normal);
			wholesumColor.R += pointColor.R * pT[i];
			wholesumColor.G += pointColor.G * pT[i];
			wholesumColor.B += pointColor.B * pT[i];
		}
	}

	double division = sqrt(m_view.nLights);
	division = 1.0 / division;
	wholesumColor.R *= division; wholesumColor.G *= division; wholesumColor.B *= division;

	// It has been intersected. Now get the color!
	// Hold on! You have to get shadow rays even for each refracted/reflected intersection!
	// if the surface is reflective, do reflection. When the ray is already inside the ball, it should not do reflection.
	Point normalizedDir = rayDir;
	normalizedDir.NormalizeVector();
	COLOR3 reflect_color(0, 0, 0), refract_color(0,0,0);
	normal.NormalizeVector();
	if(shader.Ks > 0 && (nInside == 1 || nInside == 2))
	{
		Point ReflectedRay = GetReflectRay(normal, normalizedDir);
		reflect_color = TraceRay(interP, ReflectedRay, shader.indexRefraction, ray_depth-1, nHitObj, nIgnoreID);
	}
	if(nInside == 2)
		nInside = -1;
	if(shader.Transmission > 0)
	{
		Point RefractedRay = GetRefractRay(normal, normalizedDir, _n1, shader.indexRefraction);
		// If it hit a backside, it should refract, but should not get color while transmitting/penetrating the back-surface.
		if(RefractedRay.x != SPACE_END || RefractedRay.y != SPACE_END || RefractedRay.z != SPACE_END)
		{
/*			if(Point::DotProduct(normal, normalizedDir) > 0)*/
/*			if(nInside == 1)
			{
				//wholesumColor.R /= 2; wholesumColor.G /= 2; wholesumColor.B /= 2;
				wholesumColor.R = wholesumColor.G = wholesumColor.B = 0;
				nInside = 0;
			}
			else	// When bInside == 0
				nInside = 1;*/
			refract_color = TraceRay(interP, RefractedRay, shader.indexRefraction, ray_depth-1, nHitObj, nIgnoreID, nInside);
		}
		else // No refraction if refractedRay sqrt returns nan.
			;
	}

	COLOR3 finalColor(0,0,0);
	finalColor.R = wholesumColor.R + reflect_color.R * shader.Ks + refract_color.R * shader.Transmission;
	finalColor.G = wholesumColor.G + reflect_color.G * shader.Ks + refract_color.G * shader.Transmission;
	finalColor.B = wholesumColor.B + reflect_color.B * shader.Ks + refract_color.B * shader.Transmission;
	return finalColor;
}

Point Raytracer::GetReflectRay(Point &normal, Point &rayDirection)
{
	Point pRayDir = rayDirection;
	pRayDir.NormalizeVector();
	normal.NormalizeVector();
	double c1 = Point::DotProduct(normal, pRayDir) * -2;
	c1 *= 2;
	Point R1 = pRayDir + (normal * c1); //R1 == Reflected Direction
	R1.NormalizeVector();
	return	R1;
}

Point Raytracer::GetRefractRay(Point &normal, Point &rayDirection, double _n1, double _n2)
{
/*	Point pRayDir = rayDirection;
	pRayDir.NormalizeVector();
	normal.NormalizeVector();
	double c1 = -Point::DotProduct(normal, pRayDir);
	double n1 = _n1;	// Are you sure it is 1.0?
	double n2 = _n2;
	double n = n1 / n2;
	double c2 = 1 - n*n*(1-c1*c1);
	if(c2 < 0)
		// No Refraction.
		return Point(SPACE_END, SPACE_END, SPACE_END);
	c2 = sqrt(1 - n*n*(1-c1*c1));
	Point term1 = pRayDir * n;
	Point term2 = normal * (n * c1 - c2);
	Point Rr = term1 + term2;	//Rr = Refracted Direction
	return	Rr;*/
	normal.NormalizeVector();
	Point rayDir(rayDirection);
	rayDir.NormalizeVector();
	double dir_dot_normal = Point::DotProduct(rayDir, normal);
	Point term1 = (rayDir - (normal * dir_dot_normal));
	term1 *= _n1; term1 *= (1.0/_n2);
	double square = dir_dot_normal * dir_dot_normal;
	double root_term = (1.0 - square) * _n1 * _n1;
	root_term = root_term / (_n2 * _n2);
	root_term = 1 - root_term;
	if(root_term < 0)
		return Point(SPACE_END, SPACE_END, SPACE_END);
	root_term = sqrt(root_term);
	Point term2 = normal * root_term;
	term2.x = -term2.x; term2.y = -term2.y; term2.z = -term2.z;
	Point t; t = term1 + term2;
	t.NormalizeVector();
	return	 t;
}






