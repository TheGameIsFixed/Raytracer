#include "raytracer.h"
#include <math.h>

// Mise en place du générateur de nombres aléatoires Mersenne Twister (distribution uniforme). Très bonne qualité de part ses nombreux avantages
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::mt19937 generator (seed);
std::uniform_real_distribution<double> uniform01(0.0, 1.0);

///
/// \brief Camera | Définition de notre caméra
///
CCamera Camera;

///
/// \brief CRayTracer::CRayTracer | Définition des paramètres de notre Ray Tracer
///
CRayTracer::CRayTracer()
{
        Samples = 2; // Sur-échantillonage de la scène : procédé permetant de réduire le crénelage/aliassage
        GISamples = 2; // Sur-échantillonage des ombres : procédé permetant de réduire le crénelage/aliassage
        SNSamples = 2; // Sur-échantillonage stochastique des ombres : Les lampes sont divisées en SNSamples * SNSamples

        AmbientOcclusionIntensity = 0.5f; // Intensité de l'occlusion ambiante

	Quads = NULL;
	Lights = NULL;
	Photons = NULL;

	LastQuad = NULL;
	LastLight = NULL;
	LastPhoton = NULL;

	QuadsCount = 0;
	LightsCount = 0;
	PhotonsCount = 0;

        AmbientOcclusion = true; // Occlusion ambiante : true = activé | false = désactivé
        BiDirectionnal = true; // Lancer de rayons Bi-directionnel : true = activé | false = désactivé

	PhotonsPerLight = 20; // Lancer de rayons Bi-directionnel : Nombre de rayons lancés par source
        MaxBounces = 3; // Lancer de rayons Bi-directionnel : Nombre maximal de rebonds. Permet de limiter les calculs et de fixer une condition d'arrêt si le principe de la roulette russe venait à faire rebondir les photons un grand nombre de fois ne servant à presque rien
}

///
/// \brief CRayTracer::Init | Initialisation du Ray Tracer
/// \return | true = bien initialisée - false = erreur lors de l'initialisation
///
bool CRayTracer::Init()
{
        InitScene();

	LastQuad = Quads + QuadsCount;
	LastLight = Lights + LightsCount;
	LastPhoton = Photons + PhotonsCount;

	return true;
}

///
/// \brief CRayTracer::Size | Définition de la taille de notre fenêtre
/// \param Width | largeur
/// \param Height | hauteur
///
void CRayTracer::Size(int Width, int Height)
{
	this->Width = Width;
	this->Height = Height;

	WidthMSamples = Width * Samples;
	HeightMSamples = Height * Samples;
	WidthMHeightMSamples2 = WidthMSamples * HeightMSamples;

	Camera.VPin[0] = 1.0f / (float)(WidthMSamples - 1);
	Camera.VPin[5] = 1.0f / (float)(HeightMSamples - 1);

	float aspect = (float)Width / (float)Height;
	float tany = tan(45.0f / 360.0f * (float)PI);

	Camera.Pin[0] = tany * aspect;
	Camera.Pin[5] = tany;
	Camera.Pin[10] = 0.0f;
	Camera.Pin[14] = -1.0f;

	Camera.CalculateRayMatrix();
}

///
/// \brief CRayTracer::Destroy | Destruction du Ray Tracer
///
void CRayTracer::Destroy()
{
	if(Quads != NULL)
	{
		delete [] Quads;
		Quads = NULL;
		QuadsCount = 0;
		LastQuad = NULL;
	}

	if(Lights != NULL)
	{
		delete [] Lights;
		Lights = NULL;
		LightsCount = 0;
		LastLight = NULL;
	}

	if(Photons != NULL)
	{
		delete [] Photons;
		Photons = NULL;
		PhotonsCount = 0;
		LastPhoton = NULL;
	}
}

///
/// \brief CRayTracer::Shadow | Permet de tester si un Quad fait de l'ombre à un point
/// \param Object | Quad à tester
/// \param Point | Point a tester
/// \param LightDirection | Direction du rayon de lumière
/// \param LightDistance | Distance de la lumière au point
/// \return | true = le quad fait de l'ombre au point - false = le quad ne fait pas d'ombre au point
///
bool CRayTracer::Shadow(void *Object, vec3 &Point, vec3 &LightDirection, float LightDistance)
{

	for(CQuad *Quad = Quads; Quad < LastQuad; Quad++)
	{
		if(Quad == Object) continue;

		if(Quad->Intersect(Point, LightDirection, LightDistance))
		{
			return true;
		}
	}

	return false;
}

///
/// \brief CRayTracer::LightIntensity | Calcul de l'intensité lumineuse en un point
/// \param Object | Quad
/// \param Point | Point
/// \param Normal | Normale
/// \param LightPosition | position de la lumière
/// \param Light | lumière
/// \param AO | facteur d'occultation ambiante
/// \return | intensité lumineuse du point
///
vec3 CRayTracer::LightIntensity(void *Object, vec3 &Point, vec3 &Normal, vec3 &LightPosition, CLight *Light, float AO)
{
	vec3 LightDirection = LightPosition - Point;

	float LightDistance2 = length2(LightDirection);
	float LightDistance = sqrt(LightDistance2);

	LightDirection *= 1.0f / LightDistance;

	float NdotLD = dot(Normal, LightDirection);

	if(NdotLD > 0.0f)
	{
		if(Light->Quad)
		{
			float LNdotLD = -dot(Light->Quad->N, LightDirection);

			if(LNdotLD > 0.0f)
			{
				if(Shadow(Object, Point, LightDirection, LightDistance) == false)
				{
					for(CQuad *Quad = Quads; Quad < LastQuad; Quad++)
					{
						if(Quad == Object)
						{
							vec3 V = normalize(Point-Camera.Position-Point);
							vec3 RLDLN = normalize(reflect(-LightDirection,Light->Quad->N));
							vec3 spec = Light->Specular * Quad->Ks * pow(std::max(dot(RLDLN,V),0.0f),Quad->Shininess);

							vec3 amb = Light->Ambient * AO * Quad->Ka;
							vec3 dif = Light->Diffuse * NdotLD * LNdotLD * Quad->Kd;

							return Light->Quad->Color * (amb + dif) + spec;
						}
					}
				}
			}
		}
	}
	return (Light->Quad->Color) * (Light->Ambient * AO);
}

///
/// \brief CRayTracer::PhotonIntensity | Calcul l'intensité lumineuse en un point (par les photons du lancer de rayon bi-directionnel)
/// \param Object | Quad
/// \param Point | Point
/// \param Normal | Normale
/// \param PhotonPosition | Position du photon
/// \param Photon | Photon
/// \param AO | facteur d'occultation ambiante
/// \return intensité lumineuse du photon au point d'impact
///
vec3 CRayTracer::PhotonIntensity(void *Object, vec3 &Point, vec3 &Normal, vec3 &PhotonPosition, CPhoton *Photon, float AO)
{
	vec3 PhotonDirection = PhotonPosition - Point;

	float PhotonDistance2 = length2(PhotonDirection);
	float PhotonDistance = sqrt(PhotonDistance2);

	PhotonDirection *= 1.0f / PhotonDistance;

	float NdotLD = dot(Normal, PhotonDirection);

	if(NdotLD > 0.0f)
	{
		float LNdotLD = -dot(Photon->Quad->N, PhotonDirection);

		if(LNdotLD > 0.0f)
		{
			if(Shadow(Object, Point, PhotonDirection, PhotonDistance) == false)
			{
				return Photon->totalEnergy * NdotLD * LNdotLD;
			}
		}
	}
	return vec3(0.0f);
}

///
/// \brief CRayTracer::AmbientOcclusionFactor | Calcul du facteur d'occlusion ambiante
/// \param Object | Quad
/// \param Point | Point
/// \param Normal | Normale
/// \return | facteur d'occultation ambiante
///
float CRayTracer::AmbientOcclusionFactor(void *Object, vec3 &Point, vec3 &Normal)
{
	float AO = 0.0f;

        for(int i = 0; i < GISamples; i++)
	{
		vec3 RandomRay = normalize(vec3(2 * uniform01(generator) - 1.0f, 2 * uniform01(generator) - 1.0f, 2 * uniform01(generator) - 1.0f));

		float NdotRR = dot(Normal, RandomRay);

		if(NdotRR < 0.0f)
		{
			RandomRay = -RandomRay;
			NdotRR = -NdotRR;
		}

		float Distance = INT_MAX, TestDistance;

		for(CQuad *Quad = Quads; Quad < LastQuad; Quad++)
		{
			if(Quad == Object) continue;

			if(Quad->Intersect(Point, RandomRay, Distance, TestDistance))
			{
				Distance = TestDistance;
			}
		}

		AO += NdotRR / (1.0f + Distance * Distance);
	}

        return 1.0f - AO * (1.0f / ((float)GISamples)) * AmbientOcclusionIntensity;
}

///
/// \brief CRayTracer::IlluminatePoint | Calcul de la couleur d'un point
/// \param Object
/// \param Point
/// \param Normal
/// \param Color
///
void CRayTracer::IlluminatePoint(void *Object, vec3 &Point, vec3 &Normal, vec3 &Color)
{
	float AO = 1.0f;

	if(AmbientOcclusion)
	{
		AO = AmbientOcclusionFactor(Object, Point, Normal);
	}

        if(LightsCount == 0)
        {
                float NdotCD = dot(Normal, normalize(Camera.Position - Point));

                if(NdotCD > 0.0f)
                {
                        Color *= 0.5f * (AO + NdotCD);
                }
                else
                {
                        Color *= 0.5f * AO;
                }
        }
        else
        {
                vec3 LightsIntensitiesSum;

                for(CLight *Light = Lights; Light < LastLight; Light++)
                {
                        if(Light->Quad)
                        {
                                for(int i = 0; i < GISamples; i++)
                                {
                                    for(int j = 0; j < SNSamples; j++)
                                    {
                                        for(int k = 0; k < SNSamples; k++)
                                        {
                                        float s = uniform01(generator);
                                        float t = uniform01(generator);

                                        vec3 RandomLightPosition = Light->Quad->ab * ((s+(float)j)/(float)SNSamples) + Light->Quad->ad * ((t+(float)k)/(float)SNSamples) + Light->Quad->a;

                                        LightsIntensitiesSum += LightIntensity(Object, Point, Normal, RandomLightPosition, Light, AO);
                                        }
                                    }
                                }
                        }
                }

                vec3 LightsIntensitiesSum2;
                if(BiDirectionnal == true)
                {
                        for(CPhoton *Photon = Photons; Photon < LastPhoton; Photon++)
                        {
                                LightsIntensitiesSum2 += PhotonIntensity(Object, Point, Normal, Photon->position, Photon, AO);
                        }
                }

                Color *= (LightsIntensitiesSum * (1.0f / (float)GISamples)/((float)SNSamples*(float)SNSamples) + LightsIntensitiesSum2);
        }
}

///
/// \brief CRayTracer::SendPhotonsLight | Méthode d'envoi des photons à partir d'une light
/// \param Light | lumière de départ
///
void CRayTracer::SendPhotonsLight(CLight *Light) {
	float s = uniform01(generator);
	float t = uniform01(generator);

	vec3 RandomLightPosition = Light->Quad->ab * s + Light->Quad->ad * t + Light->Quad->a;

	vec3 RandomRayNormalised = vec3(0,0,0);

	while (dot(Light->Quad->N,RandomRayNormalised)<=0)
	{
		double theta = 2 * PI * uniform01(generator);
		double phi = acos(1 - 2 * uniform01(generator));
		double x = sin(phi) * cos(theta);
		double y = sin(phi) * sin(theta);
		double z = cos(phi);
		RandomRayNormalised = normalize(vec3(x,y,z));
	}

	CData Data;
	CPhoton Photon;
	Photon.direction = RandomRayNormalised;
	Photon.energy = Light->Quad->Color / PhotonsPerLight;

	for(CQuad *Quad = Quads; Quad < LastQuad; Quad++)
	{
		if(Quad->Intersect(RandomLightPosition, RandomRayNormalised, Data.Distance, Data.TestDistance, Data.TestPoint))
		{
			Data.Point = Data.TestPoint;
			Data.Distance = Data.TestDistance;
			Data.Quad = Quad;
		}
	}



	if(Data.Quad)
	{
		Photon.position = Data.Point;
		float SolidAngle = cos(dot(Data.Quad->N,RandomRayNormalised));
		if (SolidAngle < 0)
		{
			SolidAngle = -SolidAngle;
		};
		Photon.distance = Data.Distance;
		Photon.totalEnergy = (1 - Data.Quad->Reflection) * Photon.energy * SolidAngle / ((Data.Distance) * (Data.Distance));
		Photon.energy = Data.Quad->Reflection * Photon.energy * SolidAngle;
		Photon.boucesNb ++;
		Photon.Quad = Data.Quad;

		Photons[PhotonsCount] = Photon;
		PhotonsCount++;

		double random = uniform01(generator);

		if(Data.Quad->Reflection > random && Photon.boucesNb<MaxBounces)
		{
			SendPhotonsPhoton(&Photon);
		}
	}
}


///
/// \brief CRayTracer::SendPhotonsPhoton | Méthdoe d'envoi des photons rebondis
/// \param PhotonInc | photon de départ
///
void CRayTracer::SendPhotonsPhoton(CPhoton *PhotonInc) {

	vec3 RandomRayNormalised = vec3(0,0,0);

	float AngleOuvertureDegre = PhotonInc->Quad->ReflectionAngle;
	float AngleOuvertureRadian = 2 * PI * AngleOuvertureDegre / 360;
	vec3 ReflectedRay = reflect(PhotonInc->direction, PhotonInc->Quad->N);
	RandomRayNormalised = normalize(RayCone(ReflectedRay,AngleOuvertureRadian));

	if (dot(PhotonInc->Quad->N,RandomRayNormalised)<0) {RandomRayNormalised=-RandomRayNormalised;}

	CData Data;
	CPhoton Photon;
	Photon.energy = PhotonInc->energy;
	Photon.direction = RandomRayNormalised;

	for(CQuad *Quad = Quads; Quad < LastQuad; Quad++)
	{
		if(Quad->Intersect(PhotonInc->position, RandomRayNormalised, Data.Distance, Data.TestDistance, Data.TestPoint))
		{
			Data.Point = Data.TestPoint;
			Data.Distance = Data.TestDistance;
			Data.Quad = Quad;
		}
	}

	if(Data.Quad)
	{
		Photon.position = Data.Point;
		float SolidAngle = cos(dot(Data.Quad->N,RandomRayNormalised));
		if (SolidAngle < 0)
		{
			SolidAngle = -SolidAngle;
		};
		Photon.distance = Data.Distance + PhotonInc->distance;
		Photon.totalEnergy = (1 - Data.Quad->Reflection) * Photon.energy * SolidAngle / ((Photon.distance) * (Photon.distance));
		Photon.energy = Data.Quad->Reflection * Photon.energy * SolidAngle;
		Photon.boucesNb ++;
		Photon.Quad = Data.Quad;

		Photons[PhotonsCount] = Photon;
		PhotonsCount++;

		double random = uniform01(generator);

		if(Data.Quad->Reflection > random && Photon.boucesNb<MaxBounces)
		{
			SendPhotonsPhoton(&Photon);
		}
	}
}


///
/// \brief CRayTracer::RayTrace
/// \param Origin | Origine du rayon
/// \param Ray | direction du rayon
/// \param Object
/// \return
///
vec3 CRayTracer::RayTrace(vec3 &Origin, const vec3 &Ray, void *Object)
{
	CData Data;

	for(CQuad *Quad = Quads; Quad < LastQuad; Quad++)
	{
		if(Quad == Object) continue;

		if(Quad->Intersect(Origin, Ray, Data.Distance, Data.TestDistance, Data.TestPoint))
		{
			Data.Point = Data.TestPoint;
			Data.Distance = Data.TestDistance;
			Data.Quad = Quad;
		}
	}
	for(CLight *Light = Lights; Light < LastLight; Light++)
	{
		if(Light->Quad)
		{
			if(Light->Quad->Intersect(Origin, Ray, Data.Distance, Data.TestDistance, Data.TestPoint))
			{
				Data.Point = Data.TestPoint;
				Data.Distance = Data.TestDistance;
				Data.Light = Light;
			}
		}
	}

	if(Data.Light)
	{
		Data.Color = Data.Light->Quad->Color;
	}
	else if(Data.Quad)
	{
		Data.Color = Data.Quad->Color;

		IlluminatePoint(Data.Quad, Data.Point, Data.Quad->N, Data.Color);
		
		if(Data.Quad->Reflection > 0.0f)
		{
			float AngleOuvertureDegre = Data.Quad->ReflectionAngle;
			float AngleOuvertureRadian = 2 * PI * AngleOuvertureDegre / 360;
			vec3 ReflectedRay = reflect(Ray, Data.Quad->N);
			ReflectedRay = RayCone(ReflectedRay,AngleOuvertureRadian);
                        Data.Color = mix(Data.Color, RayTrace(Data.Point, ReflectedRay, Data.Quad), Data.Quad->Reflection);
		}

	}

	return Data.Color;
}

///
/// \brief CRayTracer::InitScene | Méthode de création et d'initialisation de la scène
/// \return
///
void CRayTracer::InitScene()
{
	QuadsCount = 32; // Nombre de Quad dans la scene

	Quads = new CQuad[QuadsCount];

	LightsCount = 2; // Nombre de Lights dans la scene

	Lights = new CLight[LightsCount];

	mat4x4 Rot = RotationMatrix(-12.125f, vec3(0.0f, 1.0f, 0.0f));
	mat4x4 R = RotationMatrix(22.5f, vec3(0.0f, 1.0f, 0.0f));
	mat4x4 R2 = RotationMatrix(-15.f, vec3(1.0f, 1.0f, 1.0f));
	vec3 V = vec3(2.0f, 0.0f, 2.0f);

	//Cube
	Quads[0] =  CQuad(R * Rot * vec3(-0.5f, -2.0f,  0.5f) + V, R * Rot * vec3( 0.5f, -2.0f,  0.5f) + V, R * Rot * vec3( 0.5f, -1.0f,  0.5f) + V, R * Rot * vec3(-0.5f, -1.0f,  0.5f) + V, vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), 0.75f);
	Quads[1] =  CQuad(R * Rot * vec3( 0.5f, -2.0f, -0.5f) + V, R * Rot * vec3(-0.5f, -2.0f, -0.5f) + V, R * Rot * vec3(-0.5f, -1.0f, -0.5f) + V, R * Rot * vec3( 0.5f, -1.0f, -0.5f) + V, vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), 0.75f);
	Quads[2] =  CQuad(R * Rot * vec3( 0.5f, -2.0f,  0.5f) + V, R * Rot * vec3( 0.5f, -2.0f, -0.5f) + V, R * Rot * vec3( 0.5f, -1.0f, -0.5f) + V, R * Rot * vec3( 0.5f, -1.0f,  0.5f) + V, vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), 0.75f);
	Quads[3] =  CQuad(R * Rot * vec3(-0.5f, -2.0f, -0.5f) + V, R * Rot * vec3(-0.5f, -2.0f,  0.5f) + V, R * Rot * vec3(-0.5f, -1.0f,  0.5f) + V, R * Rot * vec3(-0.5f, -1.0f, -0.5f) + V, vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), 0.75f);
	Quads[4] =  CQuad(R * Rot * vec3(-0.5f, -1.0f,  0.5f) + V, R * Rot * vec3( 0.5f, -1.0f,  0.5f) + V, R * Rot * vec3( 0.5f, -1.0f, -0.5f) + V, R * Rot * vec3(-0.5f, -1.0f, -0.5f) + V, vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), 0.75f);
	Quads[5] =  CQuad(R * Rot * vec3(-0.5f, -2.0f, -0.5f) + V, R * Rot * vec3( 0.5f, -2.0f, -0.5f) + V, R * Rot * vec3( 0.5f, -2.0f,  0.5f) + V, R * Rot * vec3(-0.5f, -2.0f,  0.5f) + V, vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), 0.75f);

	//Cobblestone
	vec3 S = vec3(-1.25f, -3.5f, 2.25f);

	Quads[6] = CQuad(Rot * R2 * (vec3(-1.f,  1.875f,  1.f) + S), Rot * R2 * (vec3( 1.f,  1.875f,  1.f) + S), Rot * R2 * (vec3( 1.f,  1.875f, -1.f) + S), Rot * R2 * (vec3(-1.f,  1.875f, -1.f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), 0.75f);
	Quads[7] = CQuad(Rot * R2 * (vec3(-1.f,  1.875f - 0.125f,  1.f) + S), Rot * R2 * (vec3( 1.f,  1.875f - 0.125f,  1.f) + S), Rot * R2 * (vec3( 1.f,  2.0f - 0.125f,  1.f) + S), Rot * R2 * (vec3(-1.f,  2.0f - 0.125f,  1.f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), 0.75f);
	Quads[8] = CQuad(Rot * R2 * (vec3( 1.f,  1.875f - 0.125f, -1.f) + S), Rot * R2 * (vec3(-1.f,  1.875f - 0.125f, -1.f) + S), Rot * R2 * (vec3(-1.f,  2.0f - 0.125f, -1.f) + S), Rot * R2 * (vec3( 1.f,  2.0f - 0.125f, -1.f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), 0.75f);
	Quads[9] = CQuad(Rot * R2 * (vec3(-1.f,  1.875f - 0.125f, -1.f) + S), Rot * R2 * (vec3(-1.f,  1.875f - 0.125f,  1.f) + S), Rot * R2 * (vec3(-1.f,  2.0f - 0.125f,  1.f) + S), Rot * R2 * (vec3(-1.f,  2.0f - 0.125f, -1.f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), 0.75f);
	Quads[10] = CQuad(Rot * R2 * (vec3( 1.f,  1.875f - 0.125f,  1.f) + S), Rot * R2 * (vec3( 1.f,  1.875f - 0.125f, -1.f) + S), Rot * R2 * (vec3( 1.f,  2.0f - 0.125f, -1.f) + S), Rot * R2 * (vec3( 1.f,  2.0f - 0.125f,  1.f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), 0.75f);
	Quads[11] = CQuad(Rot * R2 * (vec3(-1.f,  1.875f - 0.125f, -1.f) + S), Rot * R2 * (vec3( 1.f,  1.875f - 0.125f, -1.f) + S), Rot * R2 * (vec3( 1.f,  1.875f - 0.125f,  1.f) + S), Rot * R2 * (vec3(-1.f,  1.875f - 0.125f,  1.f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), 0.75f);
		
	//Room
	//SolD
        Quads[12] =  CQuad(Rot * vec3(-0.0f, -2.0f,  5.0f), Rot * vec3( 4.0f, -2.0f,  5.0f), Rot * vec3( 4.0f, -2.0f, -4.0f), Rot * vec3(-0.0f, -2.0f, -4.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	//SolG
        Quads[13] =  CQuad(Rot * vec3(-4.0f, -2.0f,  5.0f), Rot * vec3( 0.0f, -2.0f,  5.0f), Rot * vec3( 0.0f, -2.0f,  0.0f), Rot * vec3(-4.0f, -2.0f,  0.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	//PlafondD
	Quads[14] =  CQuad(Rot * vec3( 0.0f,  2.0f, -4.0f), Rot * vec3( 4.0f,  2.0f, -4.0f), Rot * vec3( 4.0f,  2.0f,  5.0f), Rot * vec3( 0.0f,  2.0f,  5.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	//PlafondG
	Quads[15] =  CQuad(Rot * vec3(-4.0f,  2.0f,  0.0f), Rot * vec3( 0.0f,  2.0f,  0.0f), Rot * vec3( 0.0f,  2.0f,  5.0f), Rot * vec3(-4.0f,  2.0f,  5.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	//MurAD
	Quads[16] = CQuad(Rot * vec3(-0.0f, -2.0f, -4.0f), Rot * vec3( 4.0f, -2.0f, -4.0f), Rot * vec3( 4.0f,  2.0f, -4.0f), Rot * vec3(-0.0f,  2.0f, -4.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f, 0.5f, 1.f);
	//MurA
	Quads[17] = CQuad(Rot * vec3( 4.0f, -2.0f,  5.0f), Rot * vec3(-4.0f, -2.0f,  5.0f), Rot * vec3(-4.0f,  2.0f,  5.0f), Rot * vec3( 4.0f,  2.0f,  5.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	//MurD
	Quads[18] = CQuad(Rot * vec3( 4.0f, -2.0f, -4.0f), Rot * vec3( 4.0f, -2.0f,  5.0f), Rot * vec3( 4.0f,  2.0f,  5.0f), Rot * vec3( 4.0f,  2.0f, -4.0f), vec3(1.0f, 0.0f, 0.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	//MurG
	Quads[19] = CQuad(Rot * vec3(-4.0f, -2.0f,  5.0f), Rot * vec3(-4.0f, -2.0f, -0.0f), Rot * vec3(-4.0f,  2.0f, -0.0f), Rot * vec3(-4.0f,  2.0f,  5.0f), vec3(0.0f, 0.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	//MurAG
        Quads[20] = CQuad(Rot * vec3(-4.0f, -2.0f,  0.0f), Rot * vec3( 0.0f, -2.0f,  0.0f), Rot * vec3( 0.0f,  2.0f,  0.0f), Rot * vec3(-4.0f,  2.0f,  0.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f, 0.75f, 0.0f);
	//MurAM
	Quads[21] = CQuad(Rot * vec3( 0.0f, -2.0f,  0.0f), Rot * vec3( 0.0f, -2.0f, -4.0f), Rot * vec3( 0.0f,  2.0f, -4.0f), Rot * vec3( 0.0f,  2.0f,  0.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);

	//Light1
	S = vec3(-1.5f, 0.0f, 2.0f);
	Quads[22] = CQuad(Rot * (vec3(-0.5f,  1.875f,  0.5f) + S), Rot * (vec3( 0.5f,  1.875f,  0.5f) + S), Rot * (vec3( 0.5f,  1.875f, -0.5f) + S), Rot * (vec3(-0.5f,  1.875f, -0.5f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	Quads[23] = CQuad(Rot * (vec3(-0.5f,  1.875f - 0.125f,  0.5f) + S), Rot * (vec3( 0.5f,  1.875f - 0.125f,  0.5f) + S), Rot * (vec3( 0.5f,  2.0f - 0.125f,  0.5f) + S), Rot * (vec3(-0.5f,  2.0f - 0.125f,  0.5f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	Quads[24] = CQuad(Rot * (vec3( 0.5f,  1.875f - 0.125f, -0.5f) + S), Rot * (vec3(-0.5f,  1.875f - 0.125f, -0.5f) + S), Rot * (vec3(-0.5f,  2.0f - 0.125f, -0.5f) + S), Rot * (vec3( 0.5f,  2.0f - 0.125f, -0.5f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	Quads[25] = CQuad(Rot * (vec3(-0.5f,  1.875f - 0.125f, -0.5f) + S), Rot * (vec3(-0.5f,  1.875f - 0.125f,  0.5f) + S), Rot * (vec3(-0.5f,  2.0f - 0.125f,  0.5f) + S), Rot * (vec3(-0.5f,  2.0f - 0.125f, -0.5f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	Quads[26] = CQuad(Rot * (vec3( 0.5f,  1.875f - 0.125f,  0.5f) + S), Rot * (vec3( 0.5f,  1.875f - 0.125f, -0.5f) + S), Rot * (vec3( 0.5f,  2.0f - 0.125f, -0.5f) + S), Rot * (vec3( 0.5f,  2.0f - 0.125f,  0.5f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	Lights[0].Quad = new CQuad(Rot * (vec3(-0.5f,  1.875f - 0.125f, -0.5f) + S), Rot * (vec3( 0.5f,  1.875f - 0.125f, -0.5f) + S), Rot * (vec3( 0.5f,  1.875f - 0.125f,  0.5f) + S), Rot * (vec3(-0.5f,  1.875f - 0.125f,  0.5f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	Lights[0].Ambient = 0.125f;
	Lights[0].Diffuse = 0.75f;
	Lights[0].Specular = 0.125f;

	//Light2
	S = vec3(1.5f, 0.0f, 2.0f);
	Quads[27] = CQuad(Rot * (vec3(-0.5f,  1.875f,  0.5f) + S), Rot * (vec3( 0.5f,  1.875f,  0.5f) + S), Rot * (vec3( 0.5f,  1.875f, -0.5f) + S), Rot * (vec3(-0.5f,  1.875f, -0.5f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	Quads[28] = CQuad(Rot * (vec3(-0.5f,  1.875f - 0.125f,  0.5f) + S), Rot * (vec3( 0.5f,  1.875f - 0.125f,  0.5f) + S), Rot * (vec3( 0.5f,  2.0f - 0.125f,  0.5f) + S), Rot * (vec3(-0.5f,  2.0f - 0.125f,  0.5f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	Quads[29] = CQuad(Rot * (vec3( 0.5f,  1.875f - 0.125f, -0.5f) + S), Rot * (vec3(-0.5f,  1.875f - 0.125f, -0.5f) + S), Rot * (vec3(-0.5f,  2.0f - 0.125f, -0.5f) + S), Rot * (vec3( 0.5f,  2.0f - 0.125f, -0.5f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	Quads[30] = CQuad(Rot * (vec3(-0.5f,  1.875f - 0.125f, -0.5f) + S), Rot * (vec3(-0.5f,  1.875f - 0.125f,  0.5f) + S), Rot * (vec3(-0.5f,  2.0f - 0.125f,  0.5f) + S), Rot * (vec3(-0.5f,  2.0f - 0.125f, -0.5f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	Quads[31] = CQuad(Rot * (vec3( 0.5f,  1.875f - 0.125f,  0.5f) + S), Rot * (vec3( 0.5f,  1.875f - 0.125f, -0.5f) + S), Rot * (vec3( 0.5f,  2.0f - 0.125f, -0.5f) + S), Rot * (vec3( 0.5f,  2.0f - 0.125f,  0.5f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	Lights[1].Quad = new CQuad(Rot * (vec3(-0.5f,  1.875f - 0.125f, -0.5f) + S), Rot * (vec3( 0.5f,  1.875f - 0.125f, -0.5f) + S), Rot * (vec3( 0.5f,  1.875f - 0.125f,  0.5f) + S), Rot * (vec3(-0.5f,  1.875f - 0.125f,  0.5f) + S), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(1.0f, 1.0f, 1.0f), vec3(0.0f, 0.0f, 0.0f), 0.0f);
	Lights[1].Ambient = 0.125f;
	Lights[1].Diffuse = 0.75f;
	Lights[1].Specular = 0.125f;


	Photons = new CPhoton[PhotonsPerLight*MaxBounces]; //Initialisation tableau contenant le maximum possible de photons

        LastQuad = Quads + QuadsCount;
        LastLight = Lights + LightsCount;

	for(CLight *Light = Lights; Light < LastLight; Light++) //Pour chaque lampe on envoie des photons
	{
		for (int send = 0; send < PhotonsPerLight; send++)
		{
			SendPhotonsLight(Light);
		}
	}

	Camera.LookAt(vec3(0.0f), vec3(0.0f, 0.0f, 8.75f));
}

///
/// \brief CRayTracer::Draw | Méthode principale de tracé du Ray Tracer. Crée une image au format ppm
///
void CRayTracer::Draw()
{
	const int dimx = Width, dimy = Height;
	int j;
	FILE *fp = fopen("Image.ppm", "wb"); /* Image ppm générée. b - binary mode */
	(void) fprintf(fp, "P6\n%d %d\n255\n", dimx, dimy);
	for (j = dimy-1; j >= 0; j--)
	{
		{
			int Y = j * Samples;

			for(int X = 0; X < WidthMSamples; X += Samples)
			{
				static unsigned char colorOutput[3];
				vec3 SamplesSum;

				for(int yy = 0; yy < Samples; yy++)
				{
					int Yyy = Y + yy;

					for(int xx = 0; xx < Samples; xx++)
					{
						vec3 Color = RayTrace(Camera.Position, normalize(Camera.RayMatrix * vec3((float)(X + xx), (float)Yyy, 0.0f)));

						SamplesSum.r += Color.r <= 0.0f ? 0.0f : Color.r >= 1.0 ? 1.0f : Color.r;
						SamplesSum.g += Color.g <= 0.0f ? 0.0f : Color.g >= 1.0 ? 1.0f : Color.g;
						SamplesSum.b += Color.b <= 0.0f ? 0.0f : Color.b >= 1.0 ? 1.0f : Color.b;

					}
                                }

                                SamplesSum.r *= (1.0f / (float)(Samples * Samples));
                                SamplesSum.g *= (1.0f / (float)(Samples * Samples));
                                SamplesSum.b *= (1.0f / (float)(Samples * Samples));


				colorOutput[0] = SamplesSum.r*255;
				colorOutput[1] = SamplesSum.g*255;
				colorOutput[2] = SamplesSum.b*255;

				(void) fwrite(colorOutput, 1, 3, fp);

			}

		}

		// Affichage d'un compteur d'avancement
		putchar('\r');
		int pc= (dimy-j+1)*100/dimy;
		int k;
		int NBAR=30;
	     
		    printf("Progression : |");
		for (k = 0; k < NBAR; k++)
		{
			int pc_tmp = (pc * NBAR) / 100;
	       
			if (k < pc_tmp)
				putchar('-');
			else if (k == pc_tmp)
				putchar('>');
			else
				putchar(' ');
		}
		printf("| %d %%", pc);
	}
	(void) fclose(fp);
	printf("\nImage Générée\n");

}
