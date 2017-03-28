# Author: Jack Taekyu Shin (shin7@umbc.edu)

SOURCE = raytracer.cpp
IMAGE_PPM = $(IMAGE_JPG:.jpg=.ppm)
TARGET = $(SOURCE:.cpp=.out)
OBJECTS = $(SOURCE:.cpp=.o)
NFF_SOURCE = myballs.nff#balls1_transparent.nff#teapot2.nff
PPM_TARGET = $(NFF_SOURCE:.nff=.ppm)
IMAGE_JPG = $(PPM_TARGET:.ppm=.jpg)
SUPERSAMPLING = -supersampling:no
HAZE = -haze:no
THREAD_NO = 4

all : $(IMAGE_JPG)
	display $(PPM_TARGET)

$(IMAGE_JPG) : $(PPM_TARGET)
	convert $(PPM_TARGET) $(IMAGE_JPG)

$(PPM_TARGET) : $(TARGET) $(NFF_SOURCE)
	./$(TARGET) $(NFF_SOURCE) $(PPM_TARGET) $(THREAD_NO) $(SUPERSAMPLING) $(HAZE)

$(TARGET) : $(OBJECTS)
	g++ -lpthread -g -o $(TARGET) $(OBJECTS)

$(OBJECTS) : $(SOURCE)
	g++ -lpthread -g -c $(SOURCE)

# delete program, object files, and annoying xemacs ~ files
clean:
	rm -f *~ *.o *.ppm *.jpg *.out *.png

