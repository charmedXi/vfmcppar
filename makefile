# VFMCPP makefile - several versions for different compilers/architectures

# no enforced vectorisation
g++ : src/filament.h src/tangle.h src/point.h
	mkdir -p bin
	g++ -fopenmp -std=c++11 -march=native -m64 -mfpmath=sse -Ofast \
	src/main.cpp src/init.cpp src/vel.cpp src/mesh.cpp src/dummy.cpp src/field.cpp \
	src/tdt.cpp src/velnl.cpp src/mesh_adjust.cpp \
	src/ring_reconnect.cpp src/self_reconnect_line.cpp src/self_reconnect.cpp \
	src/reconnect_control.cpp src/fromfile.cpp src/loopkiller.cpp \
	-o bin/vfmcpp

gprof : src/filament.h src/tangle.h src/point.h
	mkdir -p bin
	g++ -pg -std=c++11 -march=native -m64 -mfpmath=sse -Ofast \
	src/main.cpp src/init.cpp src/vel.cpp src/mesh.cpp src/dummy.cpp src/field.cpp \
	src/tdt.cpp src/velnl.cpp src/mesh_adjust.cpp \
	src/ring_reconnect.cpp src/self_reconnect_line.cpp src/self_reconnect.cpp \
	src/reconnect_control.cpp src/fromfile.cpp src/loopkiller.cpp \
	-o bin/vfmcpp

icpc : src/filament.h src/tangle.h src/point.h
		mkdir -p bin
		/opt/intel/bin/icpc -xCORE-AVX2 -openmp -std=c++11 -O3 \
		src/main.cpp src/init.cpp src/vel.cpp src/mesh.cpp src/dummy.cpp src/field.cpp \
		src/tdt.cpp src/velnl.cpp src/mesh_adjust.cpp \
		src/ring_reconnect.cpp src/self_reconnect_line.cpp src/self_reconnect.cpp \
		src/reconnect_control.cpp src/fromfile.cpp src/loopkiller.cpp \
		-o bin/vfmcpp

