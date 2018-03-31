OBJECTS     = fluids.o
CFILES      = $(OBJECTS:.o=.c)
EXECFILE    = smoke
FFTW        = ./fftw-2.1.5
INCLUDEDIRS = -I$(FFTW)/include/
LIBDIRS     = $(FFTW)/lib
LIBS        = -framework GLUT -framework OpenGL -lrfftw -lfftw -lglui
CFLAGS      = -O2 -Wall -pipe
LINKFLAGS   = 
CC          = g++

.SILENT:

all: $(EXECFILE)

$(EXECFILE): $(OBJECTS)
		$(CC) $(LINKFLAGS) $(OBJECTS) -o $(EXECFILE) -L$(LIBDIRS) $(LIBS)

.c.o: $$@.c $$@.h
		$(CC) $(CFLAGS) $(INCLUDEDIRS) -c  $<

clean:
		-rm -rf $(OBJECTS) $(EXECFILE)

depend:
		$(CC) -MM $(CFILES) > make.dep

-include make.dep