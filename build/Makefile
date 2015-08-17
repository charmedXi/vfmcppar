.PHONY: all clean realclean

all:: vfmcpp

vpath %.cpp ../src

# optimization and language dialect flags

CXXFLAGS = -std=c++11
CXXFLAGS += -g -O3

# guess compiler from --version output and add compiler-specific
# options

compiler_key := $(shell $(CXX) --version | sed -e 's/ .*//' -e 1q)

ifeq ($(compiler_key),icpc)
CXXFLAGS += -xhost -openmp
LDFLAGS += -openmp
else
CXXFLAGS += -march=native -fopenmp
LDFLAGS += -fopenmp
endif


sources = main.cpp init.cpp vel.cpp mesh.cpp dummy.cpp field.cpp tdt.cpp velnl.cpp mesh_adjust.cpp ring_reconnect.cpp self_reconnect_line.cpp self_reconnect.cpp reconnect_control.cpp fromfile.cpp loopkiller.cpp

CPPFLAGS += -I../src

# automatic dependencies with gcc, clang, icc:
CPPFLAGS += -MMD -MP
depfiles=$(patsubst %.cpp,obj/%.d,$(sources))

-include $(depfiles)

objfiles = $(patsubst %.cpp,obj/%.o,$(sources))
obj/%.o: %.cpp
	$(CXX) -c -o $@ $(CPPFLAGS) $(CXXFLAGS) $<

vfmcpp: $(objfiles)
	$(CXX) -o $@ $(LDFLAGS) $^ $(LOADLIBES) $(LDLIBS)

clean:
	rm -f $(objfiles)

realclean: clean
	rm -f $(depfiles)
	rm -f vfmcpp

