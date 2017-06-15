#include "photon.h"

///
/// \brief CPhoton::CPhoton
///
CPhoton::CPhoton()
{
	position = vec3(0.0f);
	energy = vec3(1.0f);
	direction = vec3(1.0f);
	totalEnergy = vec3(0.0f);
	distance = 0.f;
	boucesNb = 0;

	Quad = NULL;
}
