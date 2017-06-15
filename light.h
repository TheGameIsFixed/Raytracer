#ifndef LIGHT_H
#define LIGHT_H

#include "quad.h"

class CQuad;

class CLight
{
	public:
		float Ambient, Diffuse, Specular;
		CQuad *Quad;

	public:
		CLight();
};


#endif
