#ifndef PHOTON_H
#define PHOTON_H

#include "quad.h"


class photon;

class CPhoton
{
	public:
		vec3 position, direction, energy, totalEnergy;
		float distance;
		int boucesNb;
		CQuad *Quad;

	public:
		CPhoton();
};

#endif
