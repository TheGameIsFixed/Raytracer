#include "camera.h"

///
/// \brief CCamera::CCamera
///
CCamera::CCamera()
{
	X = vec3(1.0, 0.0, 0.0);
	Y = vec3(0.0, 1.0, 0.0);
	Z = vec3(0.0, 0.0, 1.0);

	Reference = vec3(0.0, 0.0, 0.0);
	Position = vec3(0.0, 0.0, 1.0);

	Bin = BiasMatrixInverse();
}

///
/// \brief CCamera::CalculateRayMatrix
///
void CCamera::CalculateRayMatrix()
{
	Vin[0] = X.x; Vin[4] = Y.x; Vin[8] = Z.x;
	Vin[1] = X.y; Vin[5] = Y.y; Vin[9] = Z.y;
	Vin[2] = X.z; Vin[6] = Y.z; Vin[10] = Z.z;

	RayMatrix = Vin * Pin * Bin * VPin;
}

///
/// \brief CCamera::LookAt
/// \param Reference
/// \param Position
///
void CCamera::LookAt(const vec3 &Reference, const vec3 &Position)
{
	this->Reference = Reference;
	this->Position = Position;

	Z = normalize(Position - Reference);
	X = normalize(cross(vec3(0.0f, 1.0f, 0.0f), Z));
	Y = cross(Z, X);
}
