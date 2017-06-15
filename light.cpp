#include "light.h"

///
/// \brief CLight::CLight
///
CLight::CLight()
{
	Ambient = 0.0f;
	Diffuse = 1.0f;
	Specular = 0.0f;

	Quad = NULL;
}
