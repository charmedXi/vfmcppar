vfmcppar
==============
vfmcppar is a parallel C++ implementation of [vfmcpp](http://www.github.com/charmedxi/vfmcpp), a simulation code used to study the dynamics of ring/ring and ring/line scattering under the vortex filament model.

Visualisation is currently performed by a 3D matplotlib script (requires latest matplotlib version), which leaves a lot to be desired. Below are some dated animations (with old version of reconnection) to give an example visualisations only.

 ![alt text](http://giant.gfycat.com/AmbitiousPlushBetafish.gif "4 ring reconnection")

 ![alt text](http://giant.gfycat.com/ZigzagDelightfulBuzzard.gif "Highly distorted string colliding with ring")

run.sh
---------

DESCRIPTION:

Run the vfmcpp code.

USAGE:

	run.sh [OPTIONS] [FILENAME]

OPTIONS:        

	-c
		Recompile the code.
	-d
		Run the code inside the gdb debugger.
	-b 
		Run benchmarking and profiling versions.
	-h
		Show usage information.

	Flags cannot be combined.

FILENAME:

	The initial condition file can be specified; defaults to "run.in".

INITIAL CONDITIONS: 

	path 
		data folder path relative to vfmcpp root folder (will be created if doesn't exist or emptied if exists)
	
	time	
		total length of time to simulate
	
	resl	
		specify the resolution of the simulation in metres, filaments will be created, default=6.28e-8 (1um radius ring w/ 100 pts) 

	ring [radius, x, y, z, alignment]
		make a closed filament, requires 5 arguments as above, where the 5th is one of [0,1,2] for axis to align to.

	line [length, x, y, z, alignment]
		make an open filament, requires 5 arguments as above.
	
	delayed_ring [radius, x, y, z, alignment, t]
		additional argument for time t in ms to delay the addition of the ring.

	lfil [path]
		reads positions and velocities of a premade (normally distorted) line from path/pos.dat and path/vel.dat respectively.

	Eext [strength, duration, direction]
		include an external electric field, requires 3 arguments where direction is either 0, 1 or 2 for x, y or z.

	q_pt [filament, point, magnitude]
		add a charge to a filament, requires filaments to be specified eariler in file, 
		filament is the relative position of the filament in the .in file, (first filament = 0), 
		point is the index of the desired charged point on the filament, and magnitude is the size of the charge.

	#
		indicates a comment, which will be ignored by vfmcpp

EXAMPLE FILE:
	
	path data/init_example
	time 1e-3
	resl 6.28e-8
	ring 1e-6 0 0 5e-6 0 
	ring 9e-7 0 0.025e-6 0 0
	Eext 10000 1e-3 0 
	q_pt 0 50 1.6e-19
