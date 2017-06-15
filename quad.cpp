#include "quad.h"

///
/// \brief CQuad::CQuad
///
CQuad::CQuad()
{
}

///
/// \brief CQuad::CQuad | Méthode d'initialisation des valeurs d'un Quad
/// \param a | Sommet 1
/// \param b | Sommet 2
/// \param c | Sommet 3
/// \param d | Sommet 4
/// \param Color | Couleur du Quad
/// \param Ka | Coefficient de réflectivité ambiante
/// \param Kd | Coefficient de réflectivité diffuse
/// \param Ks | Coefficient de réflectivité spéculaire
/// \param Shininess | Coefficient de brillance
/// \param Reflection | Coefficient de réflection
/// \param ReflectionAngle | Angle du conne de réflexion
///
CQuad::CQuad(const vec3 &a, const vec3 &b, const vec3 &c, const vec3 &d, const vec3 &Color, const vec3 &Ka, const vec3 &Kd, const vec3 &Ks, float Shininess, float Reflection, float ReflectionAngle) : a(a), b(b), c(c), d(d), Color(Color), Ka(Ka), Kd(Kd), Ks(Ks), Shininess(Shininess), Reflection(Reflection), ReflectionAngle(ReflectionAngle)
{
	ab = b - a;
	ad = d - a;
	m = (a + b + c + d) / 4.0f;

	T = normalize(b - a);
	N = normalize(cross(b - a, c - a));

	D = -dot(N, a);

	N1 = normalize(cross(N, b - a));
	D1 = -dot(N1, a);

	N2 = normalize(cross(N, c - b));
	D2 = -dot(N2, b);

	N3 = normalize(cross(N, d - c));
	D3 = -dot(N3, c);

	N4 = normalize(cross(N, a - d));
	D4 = -dot(N4, d);
}

///
/// \brief CQuad::Inside | Méthode pour tester si un point est dans un Quad
/// \param Point | Point a tester
/// \return
///
bool CQuad::Inside(vec3 &Point)
{
	if(dot(N1, Point) + D1 < 0.0f) return false;
	if(dot(N2, Point) + D2 < 0.0f) return false;
	if(dot(N3, Point) + D3 < 0.0f) return false;
	if(dot(N4, Point) + D4 < 0.0f) return false;

	return true;
}

///
/// \brief CQuad::Intersect | Méthode pour tester si un rayon intersecte un Quad
/// \param Origin | position d'origine du rayon
/// \param Ray | direction du rayon
/// \param MaxDistance | longueur maximale à tester du rayon
/// \param Distance | longueur du rayon
/// \param Point | Point à tester
/// \return
///
bool CQuad::Intersect(vec3 &Origin, const vec3 &Ray, float MaxDistance, float &Distance, vec3 &Point)
{

	float NdotR = -dot(N, Ray);

	if(NdotR >= 0.0f)
	{

		Distance = (dot(N, Origin) + D) / NdotR;

		if(Distance >= 0.0f && Distance < MaxDistance)
		{

			Point = Ray * Distance + Origin;

			return Inside(Point);
		}
	}


	return false;
}

///
/// \brief CQuad::Intersect | Méthode pour tester si un rayon intersecte un Quad
/// \param Origin | position d'origine du rayon
/// \param Ray | direction du rayon
/// \param MaxDistance | longueur maximale à tester du rayon
/// \param Distance | longueur du rayon
/// \return
///
bool CQuad::Intersect(vec3 &Origin, const vec3 &Ray, float MaxDistance, float &Distance)
{
	float NdotR = -dot(N, Ray);

	if(NdotR > 0.0f)
	{
		Distance = (dot(N, Origin) + D) / NdotR;

		if(Distance >= 0.0f && Distance < MaxDistance)
		{
			vec3 Point;
			Point = Ray * Distance + Origin;
			return Inside(Point);
		}
	}

	return false;
}

///
/// \brief CQuad::Intersect | Méthode pour tester si un rayon intersecte un Quad
/// \param Origin | position d'origine du rayon
/// \param Ray | direction du rayon
/// \param MaxDistance | longueur maximale à tester du rayon
/// \return
///
bool CQuad::Intersect(vec3 &Origin, const vec3 &Ray, float MaxDistance)
{
	float NdotR = -dot(N, Ray);

	if(NdotR > 0.0f)
	{
		float Distance = (dot(N, Origin) + D) / NdotR;

		if(Distance >= 0.0f && Distance < MaxDistance)
		{
			vec3 Point;
			Point = Ray * Distance + Origin;
			return Inside(Point);
		}
	}

	return false;
}

