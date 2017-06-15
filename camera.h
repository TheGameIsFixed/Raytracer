#ifndef CAMERA_H
#define CAMERA_H

#include "data.h"

class CCamera
{
	public:
		vec3 X, Y, Z, Reference, Position;
		mat4x4 Vin, Pin, Bin, VPin, RayMatrix;

	public:
		CCamera();

		void CalculateRayMatrix();
		void LookAt(const vec3 &Reference, const vec3 &Position);
};

#endif
