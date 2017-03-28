# Author: Jack Taekyu Shin (shin7@umbc.edu)

SOURCE = raytracer.cpp
IMAGE_PPM = $(IMAGE_JPG:.jpg=.ppm)
TARGET = $(SOURCE:.cpp=.out)
OBJECTS = $(SOURCE:.cpp=.o)
NFF_SOURCE = balls3.nff#balls3.nff#single_triangle.nff#tetra3.nff#balls3.nff#single_triangle.nff#single_sphere.nff#balls3.nff
PPM_TARGET = $(NFF_SOURCE:.nff=.ppm)
IMAGE_JPG = $(PPM_TARGET:.ppm=.jpg)
 
all : $(IMAGE_JPG)
	display $(IMAGE_JPG)

$(IMAGE_JPG) : $(PPM_TARGET)
	convert $(PPM_TARGET) $(IMAGE_JPG)

$(PPM_TARGET) : $(TARGET) $(NFF_SOURCE)
	$(TARGET) $(NFF_SOURCE) $(PPM_TARGET)

$(TARGET) : $(OBJECTS)
	g++ -lpthread -g -o $(TARGET) $(OBJECTS)

$(OBJECTS) : $(SOURCE)
	g++ -lpthread -g -c $(SOURCE)

# delete program, object files, and annoying xemacs ~ files
clean:
	rm -f *~ *.o *.ppm *.jpg *.out

