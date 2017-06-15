#include "maths.h"

///
/// \brief mat4x4::mat4x4
///
mat4x4::mat4x4()
{
	M[0] = 1.0f; M[4] = 0.0f; M[8] = 0.0f; M[12] = 0.0f;
	M[1] = 0.0f; M[5] = 1.0f; M[9] = 0.0f; M[13] = 0.0f;
	M[2] = 0.0f; M[6] = 0.0f; M[10] = 1.0f; M[14] = 0.0f;
	M[3] = 0.0f; M[7] = 0.0f; M[11] = 0.0f; M[15] = 1.0f;
}

///
/// \brief mat4x4::~mat4x4
///
mat4x4::~mat4x4()
{
}

///
/// \brief mat4x4::mat4x4
/// \param Matrix
///
mat4x4::mat4x4(const mat4x4 &Matrix)
{
	for(int i = 0; i < 16; i++)
	{
		M[i] = Matrix.M[i];
	}
}

///
/// \brief mat4x4::operator =
/// \param Matrix
/// \return
///
mat4x4& mat4x4::operator = (const mat4x4 &Matrix)
{
	for(int i = 0; i < 16; i++)
	{
		M[i] = Matrix.M[i];
	}

	return *this;
}

///
/// \brief mat4x4::operator []
/// \param Index
/// \return
///
float& mat4x4::operator [] (int Index)
{
	return M[Index];
}

///
/// \brief mat4x4::operator &
/// \return
///
float* mat4x4::operator & ()
{
	return (float*)this;
}

///
/// \brief operator *
/// \param Matrix1
/// \param Matrix2
/// \return
///
mat4x4 operator * (const mat4x4 &Matrix1, const mat4x4 &Matrix2)
{
	mat4x4 Matrix3;

	Matrix3.M[0] = Matrix1.M[0] * Matrix2.M[0] + Matrix1.M[4] * Matrix2.M[1] + Matrix1.M[8] * Matrix2.M[2] + Matrix1.M[12] * Matrix2.M[3];
	Matrix3.M[1] = Matrix1.M[1] * Matrix2.M[0] + Matrix1.M[5] * Matrix2.M[1] + Matrix1.M[9] * Matrix2.M[2] + Matrix1.M[13] * Matrix2.M[3];
	Matrix3.M[2] = Matrix1.M[2] * Matrix2.M[0] + Matrix1.M[6] * Matrix2.M[1] + Matrix1.M[10] * Matrix2.M[2] + Matrix1.M[14] * Matrix2.M[3];
	Matrix3.M[3] = Matrix1.M[3] * Matrix2.M[0] + Matrix1.M[7] * Matrix2.M[1] + Matrix1.M[11] * Matrix2.M[2] + Matrix1.M[15] * Matrix2.M[3];

	Matrix3.M[4] = Matrix1.M[0] * Matrix2.M[4] + Matrix1.M[4] * Matrix2.M[5] + Matrix1.M[8] * Matrix2.M[6] + Matrix1.M[12] * Matrix2.M[7];
	Matrix3.M[5] = Matrix1.M[1] * Matrix2.M[4] + Matrix1.M[5] * Matrix2.M[5] + Matrix1.M[9] * Matrix2.M[6] + Matrix1.M[13] * Matrix2.M[7];
	Matrix3.M[6] = Matrix1.M[2] * Matrix2.M[4] + Matrix1.M[6] * Matrix2.M[5] + Matrix1.M[10] * Matrix2.M[6] + Matrix1.M[14] * Matrix2.M[7];
	Matrix3.M[7] = Matrix1.M[3] * Matrix2.M[4] + Matrix1.M[7] * Matrix2.M[5] + Matrix1.M[11] * Matrix2.M[6] + Matrix1.M[15] * Matrix2.M[7];

	Matrix3.M[8] = Matrix1.M[0] * Matrix2.M[8] + Matrix1.M[4] * Matrix2.M[9] + Matrix1.M[8] * Matrix2.M[10] + Matrix1.M[12] * Matrix2.M[11];
	Matrix3.M[9] = Matrix1.M[1] * Matrix2.M[8] + Matrix1.M[5] * Matrix2.M[9] + Matrix1.M[9] * Matrix2.M[10] + Matrix1.M[13] * Matrix2.M[11];
	Matrix3.M[10] = Matrix1.M[2] * Matrix2.M[8] + Matrix1.M[6] * Matrix2.M[9] + Matrix1.M[10] * Matrix2.M[10] + Matrix1.M[14] * Matrix2.M[11];
	Matrix3.M[11] = Matrix1.M[3] * Matrix2.M[8] + Matrix1.M[7] * Matrix2.M[9] + Matrix1.M[11] * Matrix2.M[10] + Matrix1.M[15] * Matrix2.M[11];

	Matrix3.M[12] = Matrix1.M[0] * Matrix2.M[12] + Matrix1.M[4] * Matrix2.M[13] + Matrix1.M[8] * Matrix2.M[14] + Matrix1.M[12] * Matrix2.M[15];
	Matrix3.M[13] = Matrix1.M[1] * Matrix2.M[12] + Matrix1.M[5] * Matrix2.M[13] + Matrix1.M[9] * Matrix2.M[14] + Matrix1.M[13] * Matrix2.M[15];
	Matrix3.M[14] = Matrix1.M[2] * Matrix2.M[12] + Matrix1.M[6] * Matrix2.M[13] + Matrix1.M[10] * Matrix2.M[14] + Matrix1.M[14] * Matrix2.M[15];
	Matrix3.M[15] = Matrix1.M[3] * Matrix2.M[12] + Matrix1.M[7] * Matrix2.M[13] + Matrix1.M[11] * Matrix2.M[14] + Matrix1.M[15] * Matrix2.M[15];

	return Matrix3;
}

///
/// \brief operator *
/// \param Matrix
/// \param Vector
/// \return
///
vec3 operator * (const mat4x4 &Matrix, const vec3 &Vector)
{
	return Matrix * vec4(Vector, 1.0f);
}

///
/// \brief operator *
/// \param Matrix
/// \param Vector
/// \return
///
vec4 operator * (const mat4x4 &Matrix, const vec4 &Vector)
{
	vec4 v;

	v.x = Matrix.M[0] * Vector.x + Matrix.M[4] * Vector.y + Matrix.M[8] * Vector.z + Matrix.M[12] * Vector.w;
	v.y = Matrix.M[1] * Vector.x + Matrix.M[5] * Vector.y + Matrix.M[9] * Vector.z + Matrix.M[13] * Vector.w;
	v.z = Matrix.M[2] * Vector.x + Matrix.M[6] * Vector.y + Matrix.M[10] * Vector.z + Matrix.M[14] * Vector.w;
	v.w = Matrix.M[3] * Vector.x + Matrix.M[7] * Vector.y + Matrix.M[11] * Vector.z + Matrix.M[15] * Vector.w;
	
	return v;
}

///
/// \brief cross
/// \param u
/// \param v
/// \return
///
vec3 cross(const vec3 &u, const vec3 &v)
{
	return vec3(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
}

///
/// \brief dot
/// \param u
/// \param v
/// \return
///
float dot(const vec3 &u, const vec3 &v)
{
	return u.x * v.x + u.y * v.y + u.z * v.z;
}

///
/// \brief length
/// \param u
/// \return
///
float length(const vec3 &u)
{
	return sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
}

///
/// \brief length2
/// \param u
/// \return
///
float length2(const vec3 &u)
{
	return u.x * u.x + u.y * u.y + u.z * u.z;
}

///
/// \brief mix
/// \param u
/// \param v
/// \param a
/// \return
///
vec3 mix(const vec3 &u, const vec3 &v, float a)
{
	return u * (1.0f - a) + v * a;
}

///
/// \brief normalize
/// \param u
/// \return
///
vec3 normalize(const vec3 &u)
{
	return u * (1.0f / sqrt(u.x * u.x + u.y * u.y + u.z * u.z));
}

///
/// \brief reflect
/// \param i
/// \param n
/// \return
///
vec3 reflect(const vec3 &i, const vec3 &n)
{
	return i - 2.0f * dot(n, i) * n;
}

///
/// \brief BiasMatrix
/// \return
///
mat4x4 BiasMatrix()
{
	mat4x4 B;

	B[0] = 0.5f; B[4] = 0.0f; B[8] = 0.0f; B[12] = 0.5f;
	B[1] = 0.0f; B[5] = 0.5f; B[9] = 0.0f; B[13] = 0.5f;
	B[2] = 0.0f; B[6] = 0.0f; B[10] = 0.5f; B[14] = 0.5f;
	B[3] = 0.0f; B[7] = 0.0f; B[11] = 0.0f; B[15] = 1.0f;

	return B;
}

///
/// \brief BiasMatrixInverse
/// \return
///
mat4x4 BiasMatrixInverse()
{
	mat4x4 BI;

	BI[0] = 2.0f; BI[4] = 0.0f; BI[8] = 0.0f; BI[12] = -1.0f;
	BI[1] = 0.0f; BI[5] = 2.0f; BI[9] = 0.0f; BI[13] = -1.0f;
	BI[2] = 0.0f; BI[6] = 0.0f; BI[10] = 2.0f; BI[14] = -1.0f;
	BI[3] = 0.0f; BI[7] = 0.0f; BI[11] = 0.0f; BI[15] = 1.0f;

	return BI;
}

///
/// \brief RotationMatrix
/// \param angle
/// \param u
/// \return
///
mat4x4 RotationMatrix(float angle, const vec3 &u)
{
	mat4x4 R;

	angle = angle / 180.0f * (float)M_PI;

	vec3 v = normalize(u);

	float c = 1.0f - cos(angle), s = sin(angle);

	R[0] = 1.0f + c * (v.x * v.x - 1.0f);
	R[1] = c * v.x * v.y + v.z * s;
	R[2] = c * v.x * v.z - v.y * s;
	R[3] = 0.0f;

	R[4] = c * v.x * v.y - v.z * s;
	R[5] = 1.0f + c * (v.y * v.y - 1.0f);
	R[6] = c * v.y * v.z + v.x * s;
	R[7] = 0.0f;

	R[8] = c * v.x * v.z + v.y * s;
	R[9] = c * v.y * v.z - v.x * s;
	R[10] = 1.0f + c * (v.z * v.z - 1.0f);
	R[11] = 0.0f;

	R[12] = 0.0f;
	R[13] = 0.0f;
	R[14] = 0.0f;
	R[15] = 1.0f;

	return R;
}


/* Cette méthode est plus réaliste que la suivante (points uniformément répartis sur la surface d'une sphère mais beaucoup plus lente) (lien: //http://corysimon.github.io/articles/uniformdistn-on-sphere/)
///
/// \brief RayCone
/// \param ReflectedRay
/// \param angleOuverture
/// \return
///
vec3 RayCone(vec3 ReflectedRay, double angleOuverture) {

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 generator (seed);
	std::uniform_real_distribution<double> uniform01(0.0, 1.0);

	ReflectedRay=normalize(ReflectedRay);

	double angle=360;
	vec3 RandomRayNormalised = vec3(0,0,0);
	while (angle > angleOuverture) {
		double theta = 2 * PI * uniform01(generator);
		double phi = acos(1 - 2 * uniform01(generator));
		double x = sin(phi) * cos(theta);
		double y = sin(phi) * sin(theta);
		double z = cos(phi);
		RandomRayNormalised = normalize(vec3(x,y,z));

		double dotRanRay = dot(ReflectedRay,RandomRayNormalised);
		angle = acos(dotRanRay)*360/PI;
	}
	ReflectedRay=normalize(RandomRayNormalised);
	return ReflectedRay;
}*/

///
/// \brief RayCone
/// \param ReflectedRay
/// \param angleOuverture
/// \return
///
vec3 RayCone(vec3 ReflectedRay, double angleOuverture) {

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 generator (seed);
	std::uniform_real_distribution<double> uniform01(0.0, 1.0);

	double RVx = uniform01(generator);
	double RVy = uniform01(generator);
	double RVz = uniform01(generator);
	vec3 RandVector = normalize(vec3(RVx,RVy,RVz));

	vec3 CrossVector = normalize(cross(ReflectedRay, RandVector));

	double s = uniform01(generator);
	double r = uniform01(generator);

	double h = cos(angleOuverture);

	double phi = 2 * PI * s;
	double z = h + ( 1 - h ) * r;
	double sinT = sqrt( 1 - z * z );
	double x = cos( phi ) * sinT;
	double y = sin( phi ) * sinT;

	vec3 perturbed = RandVector * x + CrossVector * y + ReflectedRay * z;

	ReflectedRay=normalize(perturbed);

	return ReflectedRay;
}
