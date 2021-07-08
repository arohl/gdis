# --- auto generated GDIS makefile

USE_GUI = YES
USE_GRISU = NO
include makefile.linux
include makefile.src

CFLAGS := $(CFLAGS) -DWITH_GUI 
INCS := $(INCS) `pkg-config --cflags gtk+-2.0 gthread-2.0 gtkglext-1.0 gmodule-2.0` 
LIBS := $(LIBS) `pkg-config --libs gtk+-2.0 gthread-2.0 gtkglext-1.0 gmodule-2.0` 

OBJ = $(SRC:.c=.o)
-include $(OBJ:.o=.d)

gdis: $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o ../bin/gdis $(LIBS)
%.o: %.c
	$(CC) $(CFLAGS) -c $*.c $(INCS)
	$(CC) $(CFLAGS) -MM $*.c $(INCS) > $*.d
clean:
	 rm -rf gdis *.o *.d
