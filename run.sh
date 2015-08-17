#!/bin/zsh

red='\e[1;37m'
blu='\e[1;34m'
yel='\e[1;33m'
gre='\e[0;36m'


usage() {
cat <<EOF

#########################################################################
#                                                                       #
#                                                                       #
#            ██╗   ██╗███████╗███╗   ███╗ ██████╗██████╗ ██████╗        #
#            ██║   ██║██╔════╝████╗ ████║██╔════╝██╔══██╗██╔══██╗       #
#            ██║   ██║█████╗  ██╔████╔██║██║     ██████╔╝██████╔╝       #
#            ╚██╗ ██╔╝██╔══╝  ██║╚██╔╝██║██║     ██╔═══╝ ██╔═══╝        #
#             ╚████╔╝ ██║     ██║ ╚═╝ ██║╚██████╗██║     ██║            #
#              ╚═══╝  ╚═╝     ╚═╝     ╚═╝ ╚═════╝╚═╝     ╚═╝            #
#                                                                       #
#               * * *            vfmcpp            * * *                #
#                                                                       #
#                         vortex filament code                          #
#                      Rory Brown & Matthew Evans                       #
#                                                                       #
#########################################################################

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


EOF
}

CPILE=0
DEBUG=0
BENCH=0
no=$#

cwd=$(pwd)

while true
do
	case "$4" in
		-h|--help) usage; shift;;
		-c|--cpile) CPILE=1; shift;;
		-d|--debug) DEBUG=1; shift;;
        -b|--bench) BENCH=1; shift;;
		-*) echo "invalid flag, showing help"; usage; shift;; 
		*) in=$4; break;;	
	esac
done

while true
do
	case "$3" in
		-h|--help) usage; shift;;
		-c|--cpile) CPILE=1; shift;;
		-d|--debug) DEBUG=1; shift;;
        -b|--bench) BENCH=1; shift;;
		-*) echo "invalid flag, showing help"; usage; shift;; 
		*) in=$3; break;;	
	esac
done

while true
do
	case "$2" in
		-h|--help) usage; shift;;
		-c|--cpile) CPILE=1; shift;;
		-d|--debug) DEBUG=1; shift;;
        -b|--bench) BENCH=1; shift;;
		-*) echo "invalid flag, showing help"; usage; shift;; 
		*) in=$2; break;;	
	esac
done

while true
do
	case "$1" in
		-h|--help) usage; shift;;
		-c|--cpile) CPILE=1; shift;;
		-d|--debug) DEBUG=1; shift;;
        -b|--bench) BENCH=1; shift;;
		-*) echo "invalid flag, showing help"; usage; shift;; 
		*) in=$1; break;;	
	esac
done

if [ ! -f "$in" ]; then
	echo "given input file does not exist, exiting..."
	exit 1
fi


dir=$(grep -v '^#' $in | grep 'path')
dir=${dir:5}
echo $dir
	
if [ -d "$dir" ]; then
	echo "data directory found, cleaning up...\n\n"
	cd "$dir"
	rm *.dat
	if [ ! -d "snapshot" ]; then
		mkdir "snapshot"
	fi
	cd $cwd
fi

if [ ! -d "$dir" ]; then
	echo "data directory not found, creating...\n\n"
	mkdir -p "$dir"
	mkdir "$dir/snapshot"
fi

if [ $CPILE -eq 1 ]; then
	echo " ${blu} -c flag specified, recompiling source...${red}"
	[ -d bin ] || mkdir bin
        [ -e bin/vfmcpp ] || ln -s build/vfmcpp bin/vfmcpp
	(cd build && make)
	echo " ${blu} success!\n\n"
fi

if [ $BENCH -eq 1 ]; then
	echo " ${blu} -b flag specified, running benchmark...${red}"
	echo " ${blu} running executable in bin/chmarking/icpc"
	cd "bin/chmarking/icpc-autovec"
fi

if [ $BENCH -eq 0 ]; then
	cd "bin"
fi

echo "  ${gre}#########################################################################"
echo "  ${gre}#                                                                       #"
echo "  ${gre}#                                                                       #"
echo "  ${gre}#         ${blu} ██╗   ██╗███████╗███╗   ███╗ ██████╗██████╗ ██████╗          ${gre}#"
echo "  ${gre}#         ${blu} ██║   ██║██╔════╝████╗ ████║██╔════╝██╔══██╗██╔══██╗         ${gre}#"
echo "  ${gre}#         ${blu} ██║   ██║█████╗  ██╔████╔██║██║     ██████╔╝██████╔╝         ${gre}#"
echo "  ${gre}#         ${blu} ╚██╗ ██╔╝██╔══╝  ██║╚██╔╝██║██║     ██╔═══╝ ██╔═══╝          ${gre}#"
echo "  ${gre}#         ${blu}  ╚████╔╝ ██║     ██║ ╚═╝ ██║╚██████╗██║     ██║              ${gre}#"
echo "  ${gre}#         ${blu}   ╚═══╝  ╚═╝     ╚═╝     ╚═╝ ╚═════╝╚═╝     ╚═╝              ${gre}#"
echo "  ${gre}#                                                                       #"
echo "  ${gre}#               ${yel}* * *            vfmcpp            * * *         ${gre}       #"
echo "  ${gre}#                                                                       ${gre}#"
echo "  ${gre}#                 ${yel}        vortex filament code                          ${gre}#"
echo "  ${gre}#                 ${yel}     Rory Brown & Matthew Evans                       ${gre}#"
echo "  ${gre}#                                                                       ${gre}#"
echo "  ${gre}#########################################################################\n${red}"


if [ $DEBUG -eq 1 ]; then
	if [ $BENCH -eq 1 ]; then
		gdb -ex run --args vfmcpp "../../../$in"
	fi
	if [ $BENCH -eq 0 ]; then
		gdb -ex run --args vfmcpp "../$in"
	fi

fi

if [ ! $DEBUG -eq 1 ]; then
	if [ $BENCH -eq 1 ]; then
		./vfmcpp "../../../$in"
	fi
	if [ $BENCH -eq 0 ]; then
		./vfmcpp "../$in"
	fi
fi
