# --- modules Makefile

# --- edit as required

CC = gcc

# shared object (linux)
SO = -shared
# shared object (OSX)
#SO = -bundle -flat_namespace -undefined suppress

CFLAGS = -fPIC -Wall
INCS = `pkg-config --cflags gtk+-2.0 gthread-2.0 gtkglext-1.0 gmodule-2.0`
LIBS = `pkg-config --libs gtk+-2.0 gthread-2.0 gtkglext-1.0 gmodule-2.0`

# --- main target
SRC1 = phonon.c
#SRC2 = dock.c

OBJ1 = $(SRC1:.c=.o)
OBJ2 = $(SRC2:.c=.o)

.c.o: 
	$(CC) $(CFLAGS) -c $< $(INCS)

# --- extra targets

all: phonon 

phonon: $(OBJ1)
	$(CC) $(SO) -o libphonon.so $(OBJ1)

dock: $(OBJ2)
	$(CC) $(SO) -o libdock.so $(OBJ2)

.PHONY: clean
clean:
	/bin/rm -f $(OBJ) lib*.so
