vfmcpp : src/filament.h src/tangle.h src/point.h
	g++ -std=c++11 -g -O3 src/main.cpp src/init.cpp src/vel.cpp src/mesh.cpp src/dummy.cpp src/field.cpp src/tdt.cpp src/velnl.cpp src/mesh_adjust.cpp src/reconnect.cpp src/fromfile.cpp src/loopkiller.cpp -o bin/vfmcpp