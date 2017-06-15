EXE_NAME = raytrace 
OBJS = main.o raytracer.o quad.o light.o camera.o data.o maths.o photon.o
CXX = g++
CXXFLAGS = -g -Wall -Wno-reorder -std=c++11
LIBS =

.c.o:
	$(CXX) $< -c $(CXXFLAGS) $(INC)

$(EXE_NAME): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(CXXFLAGS) $(LIBS)
clean:
	rm -f $(OBJS) $(EXE_NAME)
