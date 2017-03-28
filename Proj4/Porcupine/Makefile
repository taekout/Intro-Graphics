# GNU make magic to get the right platform-dependent flags for compiling. 
# For example, on linux this will include Defs.Linux, on mac it will include
# Defs.Darwin and on window with mingw it will include Defs.MINGW32_NT-6.1
include Defs.$(shell uname)

# files and intermediate files we create
PROG=porcupine
SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)

# rules for building -- ordered from final output to original .c for no
# particular reason other than that the first rule is the default

all : $(PROG)
	$(PROG)

# program from .o files
$(PROG): $(OBJS)
	$(CXX) -g3 -O0 -o $(PROG) $(OBJS) $(LDFLAGS) $(LDLIBS)

# .o from .cpp
%.o: %.cpp
	$(CXX) -g3 -O0 -c -o $@ $< $(CXXFLAGS)

# remove everything but the program
clean:
	rm -f *~ *.o

# remove everything including program
clobber: clean
	rm -f $(PROG)$(EXE)
