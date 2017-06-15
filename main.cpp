#include "raytracer.h"

CRayTracer RayTracer;

///
/// \brief main
/// \param argc
/// \param argv
/// \return
///
int main(int argc, char** argv)
{
	RayTracer.Size(800,600);
	RayTracer.Init();
	RayTracer.Draw();
	RayTracer.Destroy();

	return 0;
}
