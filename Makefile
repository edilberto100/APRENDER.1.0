CC = gcc
INCLUDE = -I/usr/local/hdf5/include -I/usr/local/cfitsio/include
CFLAGS = $(INCLUDE)
LDFLAGS = -lm -lhdf5 -lcfitsio
LIBS = -L/usr/local/cfitsio/lib/ -L/usr/local/hdf5/lib

APRENDER_H = aprender_fits.c aprender_reduction.c aprender_hdf5.c aprender_search.c aprender_photometry.c aprender_field.c

APRENDER = aprender.c $(APRENDER_H) 

OBJS_Field = $(APRENDER:.c=.o)
DEPS_Field = $(APRENDER:.c=.d)

TARGET_Field = APRENDER.6.0

all: 
	$(CC) $(APRENDER) $(CFLAGS) $(LDFLAGS) $(LIBS) -o $(TARGET_Field)

clean:
	$(RM) $(OBJS_Field) $(OBJS_Star) core *~ *.bak

empty: clean
	$(RM) $(DEPS)

.PHONY: all install clean empty

