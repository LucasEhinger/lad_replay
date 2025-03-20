# Generated automatically from Makefile.in by configure.
# This Makefile is used to create several test cases for Tcl's load
# command.  It also illustrates how to take advantage of configuration
# exported by Tcl to set up Makefiles for shared libraries.
# RCS: @(#) $Id: Makefile.in,v 1.6 1999/09/21 06:37:32 hobbs Exp $

TCL_DBGX =		
CC = gcc -pipe
#SCALERLIBS = -L../../coda25/crl/RPC -lscaser
SCALERLIBS = -L.  -lscaser
#LIBS =			-L/usr/lib64 -ltclstub8.6${TCL_DBGX} -ldl  -lieee $(SCALERLIBS) -lm -lc
LIBS =			-L/usr/lib64 -ltclstub8.6${TCL_DBGX} -ldl $(SCALERLIBS) -lm -lc
AC_FLAGS =		
SHLIB_CFLAGS =		-fPIC
SHLIB_LD =		${CC} -shared -Wl,-soname,${@}
SHLIB_SUFFIX =		.so
SHLIB_VERSION =		
SRC_DIR =		.
TCL_VERSION=		8.3

CFLAGS = -g
#CC_SWITCHES = $(CFLAGS) -I${SRC_DIR}/../../generic -DTCL_MEM_DEBUG \
#	${SHLIB_CFLAGS} -DUSE_TCL_STUBS ${AC_FLAGS} -I./include \
#	-I$(TCLINCLUDE_DIR)
CC_SWITCHES = $(CFLAGS) -I${SRC_DIR}/../../generic -DTCL_MEM_DEBUG \
	${SHLIB_CFLAGS} -DUSE_TCL_STUBS ${AC_FLAGS} -I./include  -I /usr/include/tirpc

all: scaler${SHLIB_SUFFIX}

pkga${SHLIB_SUFFIX}: $(SRC_DIR)/pkga.c
	$(CC) -c $(CC_SWITCHES) $(SRC_DIR)/pkga.c
	${SHLIB_LD} -o pkga${SHLIB_SUFFIX} pkga.o ${LIBS}

pkgb${SHLIB_SUFFIX}: $(SRC_DIR)/pkgb.c
	$(CC) -c $(CC_SWITCHES) $(SRC_DIR)/pkgb.c
	${SHLIB_LD} -o pkgb${SHLIB_SUFFIX} pkgb.o ${LIBS}

scaler${SHLIB_SUFFIX}: $(SRC_DIR)/scaler.c
	$(CC) -c $(CC_SWITCHES) $(SRC_DIR)/scaler.c
	${SHLIB_LD} -o scaler${SHLIB_SUFFIX} scaler.o ${LIBS}

genscaler${SHLIB_SUFFIX}: $(SRC_DIR)/genscaler.c $(SRC_DIR)/genscalib.o
	$(CC) -c $(CC_SWITCHES) $(SRC_DIR)/genscaler.c
	${SHLIB_LD} -o genscaler${SHLIB_SUFFIX} genscaler.o $(SRC_DIR)/genscalib.o ${LIBS}

clean:
	rm -f *.o *${SHLIB_SUFFIX} config.cache config.log config.status lib.exp

distclean: clean
	rm -f Makefile
