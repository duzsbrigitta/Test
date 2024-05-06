# -----------------------------------------------------------------


SHELL = /bin/bash

includedir   = /usr/local/include
libdir       = /usr/local/lib

LIBTOOL      = $(SHELL) /uni-mainz.de/homes/bridzs/Documents/sundials-2.5.0/libtool
LIBTOOL_DEPS = config/ltmain.sh
CPP      = cc -E
CC       = cc
CFLAGS   = -g -m64 -Ofast -flto 
LIBS     = -lm 

TARGET = general-pc

SUNDIALS_INCS = -I$/uni-mainz.de/homes/bridzs/Documents/sundials-2.5.0/include

SUNDIALS_LIBS = /uni-mainz.de/homes/bridzs/Documents/sundials-2.5.0/src/cvode/libsundials_cvode.la            \
	        /uni-mainz.de/homes/bridzs/Documents/sundials-2.5.0/src/nvec_ser/libsundials_nvecserial.la

all:
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(SUNDIALS_INCS) $(CFLAGS) -c $(TARGET).c -o $(TARGET).o ; \
	$(LIBTOOL) --mode=link $(CC) -o $(TARGET) $(TARGET).o $(CFLAGS) $(SUNDIALS_LIBS) $(LIBS); \


 	
