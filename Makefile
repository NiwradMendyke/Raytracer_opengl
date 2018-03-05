# Makefile
# CSCI 480
# Assignment 3

# we assume the pic directory locates one level above,
# change PIC_PATH if this is not the case
PIC_PATH = $(abspath $(CURDIR)/../pic)

INCLUDE = -I$(PIC_PATH)
LIBRARIES = -L$(PIC_PATH) -framework OpenGL -framework GLUT -lpicio $(PIC_PATH)/libjpeg.a

COMPILER = g++ -std=c++11
COMPILERFLAGS = -O3 $(INCLUDE) -Wno-deprecated-declarations

PROGRAM = assign3
SOURCE = assign3.cpp
OBJECT = assign3.o

.cpp.o:
	$(COMPILER) -c $(COMPILERFLAGS) $<

all: $(PROGRAM)

$(PROGRAM): $(OBJECT)
	$(COMPILER) $(COMPILERFLAGS) -o $(PROGRAM) $(OBJECT) $(LIBRARIES)

clean:
	-rm -rf core *.o *~ "#"*"#" $(PROGRAM)
