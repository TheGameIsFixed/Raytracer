#ifndef MATHS_H
#define MATHS_H

#include <chrono>
#include <random>
# define PI           (atan(1) * 4)  /* pi */

class vec3
{
	public:
		union{
			struct{float x, y, z;};
			struct{float s, t, p;};
			struct{float r, g, b;};
		};
		vec3() : x(0.0f), y(0.0f), z(0.0f){}
		vec3(float num) : x(num), y(num), z(num){}
		vec3(float x, float y, float z) : x(x), y(y), z(z){}
		vec3& operator = (const vec3 &u){x = u.x; y = u.y; z = u.z; return 	*this;}
		vec3 operator - (){return vec3(-x, -y, -z);}
		vec3& operator += (float num){x += num; y += num; z += num; return *this;}
		vec3& operator += (const vec3 &u){x += u.x; y += u.y; z += u.z; return *this;}
		vec3& operator -= (float num){x -= num; y -= num; z -= num; return *this;}
		vec3& operator -= (const vec3 &u){x -= u.x; y -= u.y; z -= u.z; return *this;}
		vec3& operator *= (float num){x *= num; y *= num; z *= num; return *this;}
		vec3& operator *= (const vec3 &u){x *= u.x; y *= u.y; z *= u.z; return *this;}
		vec3& operator /= (float num){x /= num; y /= num; z /= num; return *this;}
		vec3& operator /= (const vec3 &u){x /= u.x; y /= u.y; z /= u.z; return *this;}
		friend vec3 operator + (const vec3 &u, float num){return vec3(u.x + num, u.y + num, u.z + num);}
		friend vec3 operator + (float num, const vec3 &u){return vec3(num + u.x, num + u.y, num + u.z);}
		friend vec3 operator + (const vec3 &u, const vec3 &v){return vec3(u.x + v.x, u.y + v.y, u.z + v.z);}
		friend vec3 operator - (const vec3 &u, float num){return vec3(u.x - num, u.y - num, u.z - num);}
		friend vec3 operator - (float num, const vec3 &u){return vec3(num - u.x, num - u.y, num - u.z);}
		friend vec3 operator - (const vec3 &u, const vec3 &v){return vec3(u.x - v.x, u.y - v.y, u.z - v.z);}
		friend vec3 operator * (const vec3 &u, float num){return vec3(u.x * num, u.y * num, u.z * num);}
		friend vec3 operator * (float num, const vec3 &u){return vec3(num * u.x, num * u.y, num * u.z);}
		friend vec3 operator * (const vec3 &u, const vec3 &v){return vec3(u.x * v.x, u.y * v.y, u.z * v.z);}
		friend vec3 operator / (const vec3 &u, float num){return vec3(u.x / num, u.y / num, u.z / num);}
		friend vec3 operator / (float num, const vec3 &u){return vec3(num / u.x, num / u.y, num / u.z);}
		friend vec3 operator / (const vec3 &u, const vec3 &v){return vec3(u.x / v.x, u.y / v.y, u.z / v.z);}
};


class vec4
{
	public:
		union{
			struct{float x, y, z, w;};
			struct{float s, t, p, q;};
			struct{float r, g, b, a;};
		};
		vec4() : x(0.0f), y(0.0f), z(0.0f), w(0.0f){}
		~vec4(){}
		vec4(float num) : x(num), y(num), z(num), w(num){}
		vec4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w){}
		vec4(const vec3 &u, float w) : x(u.x), y(u.y), z(u.z), w(w){}
		vec4(const vec4 &u) : x(u.x), y(u.y), z(u.z), w(u.w){}
		vec4& operator = (const vec4 &u){x = u.x; y = u.y; z = u.z; w = u.w; return *this;}
		vec4 operator - (){return vec4(-x, -y, -z, -w);}
		float* operator & (){return (float*)this;}
		operator vec3 (){return *(vec3*)this;}
		vec4& operator += (float num){x += num; y += num; z += num; w += num; return *this;}
		vec4& operator += (const vec4 &u){x += u.x; y += u.y; z += u.z; w += u.w; return *this;}
		vec4& operator -= (float num){x -= num; y -= num; z -= num; w -= num; return *this;}
		vec4& operator -= (const vec4 &u){x -= u.x; y -= u.y; z -= u.z; w -= u.w; return *this;}
		vec4& operator *= (float num){x *= num; y *= num; z *= num; w *= num; return *this;}
		vec4& operator *= (const vec4 &u){x *= u.x; y *= u.y; z *= u.z; w *= u.w; return *this;}
		vec4& operator /= (float num){x /= num; y /= num; z /= num; w /= num; return *this;}
		vec4& operator /= (const vec4 &u){x /= u.x; y /= u.y; z /= u.z; w /= u.w; return *this;}
		friend vec4 operator + (const vec4 &u, float num){return vec4(u.x + num, u.y + num, u.z + num, u.w + num);}
		friend vec4 operator + (float num, const vec4 &u){return vec4(num + u.x, num + u.y, num + u.z, num + u.w);}
		friend vec4 operator + (const vec4 &u, const vec4 &v){return vec4(u.x + v.x, u.y + v.y, u.z + v.z, u.w + v.w);}
		friend vec4 operator - (const vec4 &u, float num){return vec4(u.x - num, u.y - num, u.z - num, u.w - num);}
		friend vec4 operator - (float num, const vec4 &u){return vec4(num - u.x, num - u.y, num - u.z, num - u.w);}
		friend vec4 operator - (const vec4 &u, const vec4 &v){return vec4(u.x - v.x, u.y - v.y, u.z - v.z, u.w - v.w);}
		friend vec4 operator * (const vec4 &u, float num){return vec4(u.x * num, u.y * num, u.z * num, u.w * num);}
		friend vec4 operator * (float num, const vec4 &u){return vec4(num * u.x, num * u.y, num * u.z, num * u.w);}
		friend vec4 operator * (const vec4 &u, const vec4 &v){return vec4(u.x * v.x, u.y * v.y, u.z * v.z, u.w * v.w);}
		friend vec4 operator / (const vec4 &u, float num){return vec4(u.x / num, u.y / num, u.z / num, u.w / num);}
		friend vec4 operator / (float num, const vec4 &u){return vec4(num / u.x, num / u.y, num / u.z, num / u.w);}
		friend vec4 operator / (const vec4 &u, const vec4 &v){return vec4(u.x / v.x, u.y / v.y, u.z / v.z, u.w / v.w);}
};


class mat4x4
{
	public:
		float M[16];
		mat4x4();
		~mat4x4();
		mat4x4(const mat4x4 &Matrix);
		mat4x4& operator = (const mat4x4 &Matrix);
		float& operator [] (int Index);
		float* operator & ();
		friend mat4x4 operator * (const mat4x4 &Matrix1, const mat4x4 &Matrix2);
		friend vec3 operator * (const mat4x4 &Matrix, const vec3 &Vector);
		friend vec4 operator * (const mat4x4 &Matrix, const vec4 &Vector);
};


vec3 cross(const vec3 &u, const vec3 &v);
float dot(const vec3 &u, const vec3 &v);
float length(const vec3 &u);
float length2(const vec3 &u);
vec3 mix(const vec3 &u, const vec3 &v, float a);
vec3 normalize(const vec3 &u);
vec3 reflect(const vec3 &i, const vec3 &n);
vec3 RayCone(vec3 ReflectedRay, double angleOuverture);

mat4x4 BiasMatrix();
mat4x4 BiasMatrixInverse();
mat4x4 RotationMatrix(float angle, const vec3 &u);

#endif
