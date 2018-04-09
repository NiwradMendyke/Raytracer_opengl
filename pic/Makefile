CC = gcc

OBJS = pic.o xpic.o ppm.o adaptcm.o jpeg.o

LIB = libpicio.a

CPPFLAGS = -I/opt/X11/include 
LDFLAGS = -L/opt/X11/lib

$(LIB): $(OBJS)
	ar cr $(LIB) $(OBJS)

clean:
	/bin/rm -f $(LIB) $(OBJS)
