#ifndef RAYTRACER_H
#define RAYTRACER_H

#include "data.h"
#include "camera.h"
#include "light.h"
#include "photon.h"

class CRayTracer
{
	private:
		int Width, LineWidth, Height, Samples, GISamples, SNSamples, WidthMSamples, HeightMSamples, WidthMHeightMSamples2;
		float AmbientOcclusionIntensity;

	protected:
		CQuad *Quads, *LastQuad;
		CLight *Lights, *LastLight;
		CPhoton *Photons, *LastPhoton;
		int QuadsCount, LightsCount, PhotonsCount, PhotonsPerLight, MaxBounces;

	public:
		bool AmbientOcclusion, BiDirectionnal;

	public:
		CRayTracer();

		bool Init();
		void Size(int Width, int Height);
		void Draw();
		void Destroy();
		vec3 RayTrace(vec3 &Origin, const vec3 &Ray, void *Object = NULL);


	protected:
		void InitScene();

	private:
		bool Shadow(void *Object, vec3 &Point, vec3 &LightDirection, float LightDistance);
		vec3 LightIntensity(void *Object, vec3 &Point, vec3 &Normal, vec3 &LightPosition, CLight *Light, float AO);
		vec3 PhotonIntensity(void *Object, vec3 &Point, vec3 &Normal, vec3 &PhotonPosition, CPhoton *Photon, float AO);
		float AmbientOcclusionFactor(void *Object, vec3 &Point, vec3 &Normal);
		void IlluminatePoint(void *Object, vec3 &Point, vec3 &Normal, vec3 &Color);
		void SendPhotonsLight(CLight *Light);
		void SendPhotonsPhoton(CPhoton *PhotonInc);

};

#endif
