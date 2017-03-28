#include <iostream>
#include <math.h>

using namespace std;


class Point{
public:
	double x;
	double y;
	double z;

public:
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
	Point & operator*=(const double t)
	{
		x *= t;
		y *= t;
		z *= t;
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

	float* GetPoint(void)
	{
//		Point ret(x, y, z);
//		return	ret;
		float *p = new float[3];
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

	static void CrossProduct(Point p1, Point p2, Point &p3)
	{
		float *M1 = p1.GetPoint();
		float *M2 = p2.GetPoint();
		p3.x = M1[1] * M2[2] - M1[2] * M2[1];
		p3.y = M1[2] * M2[0] - M1[0] * M2[2];
		p3.z = M1[0] * M2[1] - M1[1] * M2[0];
	}

	static double DotProduct(const Point p1, const Point p2)
	{
		return	p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
	}
	void Print()
	{
		cout << x << " " << y << " " << z << endl;
	}
};

void testfunc(Point p)
{
	p . x = 100;
	p . y = 100;
	p . z = 100;
}

void func(Point &tmp)
{
	tmp.x += 100;
}

int main(void)
{
	Point a(2, 1.3, 1.7);
	Point c;
	Point::CrossProduct(a, Point(0,1,0), c);
	c.NormalizeVector();
	Point d;
	Point::CrossProduct(a, c, d);
	d.NormalizeVector();
	a = a * -1;
	Point e[4];
	Point f[4];
	e[0] = (d * -1.5) + (c*0.3);
	e[1] = (d * -1.5) + (c*-0.3);
	e[2] = (d * 1.5) + (c*-0.3);
	e[3] = (d * 1.5) + (c*0.3);
	f[0] = (d * 0.3) + (c * 1.5);
	f[1] = (d * 0.3) + (c * -1.5);
	f[2] = (d * -0.3) + (c * 1.5);
	f[3] = (d * -0.3) + (c * 1.5);
	for(int i = 0 ; i < 4 ; i++)
	{
		e[i] = e[i] + a;
		f[i] = f[i] + a;
	}
	for(int i = 0 ; i < 4 ; i++)
	{
		e[i].Print();
		f[i].Print();
	}
	return	0;
}