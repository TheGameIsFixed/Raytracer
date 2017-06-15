#ifndef UTILS_H
#define UTILS_H

#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <iostream>
#include "maths.h"
#include "quad.h"
#include "light.h"

class CQuad;
class CLight;

class CData
{
public:
	float Distance, TestDistance;
	vec3 Color, Point, TestPoint;
	CQuad *Quad;
	CLight *Light;

public:
	CData();
};

#endif
