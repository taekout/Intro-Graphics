# path to GLUT and GL
GLUTROOT = /afs/umbc.edu/users/o/l/olano/home/public/freeglut
GLROOT = /usr/local
LIBGLU = /usr/local/lib/libGLU.so
LIBGL = /usr/local/lib/libGL.so

NFF_FILE = level.nff

# flags for C compiler
#   glut under $(GLROOT)
#   -g to enable debugging info
CFLAGS = -I$(GLUTROOT)/include -I$(GLROOT)/include -g

# flags for C++ compiler -- same as C compiler so C++ can find includes
CXXFLAGS = $(CFLAGS)

# where to find libraries
#   -L to find when linking and 
#   -rpath to find shared libraries at run-time
LDFLAGS = -L$(GLUTROOT)/lib -Wl,-rpath,$(GLUTROOT)/lib,-rpath,$(GLROOT)/lib

# names of libraries to use:
#   glut for cross-platform window creation and interaction
#   GLU for a set of utilities for GL (some of which are used by GLUT)
#   GL for OpenGL itself
LDLIBS = -lglut $(LIBGL) $(LIBGLU) -lm

# files and intermediate files we create
OBJS  = example.o view.o motion.o draw.o
PROG  = example

# rules for building -- ordered from final output to original .c for no
# particular reason other than that the first rule is the default

# program from .o files
all : $(PROG)
	$(PROG) $(NFF_FILE)
$(PROG): $(OBJS)
	$(CXX) -o $(PROG) $(OBJS) $(LDFLAGS) $(LDLIBS)

# any .o from .c uses built-in rule


# remove everything but the program
clean:
	rm -f *~ *.o

# remove everything including program
clobber: clean
	rm -f $(PROG)
