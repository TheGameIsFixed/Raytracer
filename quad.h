#ifndef QUAD_H
#define QUAD_H

#include "data.h"

class CQuad
{
	public:
		float Reflection, ReflectionAngle, D, D1, D2, D3, D4;
		vec3 a, b, c, d, Color, ab, ad, m, T, N, N1, N2, N3, N4;
		vec3 Ka, // Coefficient de réflectivité ambiante
		      Kd, // Coefficient de réflectivité diffuse
		      Ks; // Coefficient de réflectivité spéculaire
		float Shininess; // Coefficient de brillance

	public:
		CQuad();
		CQuad(const vec3 &a, const vec3 &b, const vec3 &c, const vec3 &d, const vec3 &Color, const vec3 &Ka, const vec3 &Kd, const vec3 &Ks, float Shininess, float Reflection = 0.0f, float ReflectionAngle = 0.0f);

		bool Inside(vec3 &Point);
		bool Intersect(vec3 &Origin, const vec3 &Ray, float MaxDistance, float &Distance, vec3 &Point);
		bool Intersect(vec3 &Origin, const vec3 &Ray, float MaxDistance, float &Distance);
		bool Intersect(vec3 &Origin, const vec3 &Ray, float MaxDistance);
};

#endif
