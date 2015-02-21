# use this after installing mpir http://www.mpir.org/ to see if you get better performance than using gmp
#INT_LIB=mpir
CXX = g++
CC = cc
CXXFLAGS = -Wall -Wextra -std=c++0x -O3 -fomit-frame-pointer  -g 

CFLAGS = -Wall -Wextra -O3 -fomit-frame-pointer  -g -std=gnu99 


OSVERSION := $(shell uname -s)
LIBS = -lcrypto -lssl -pthread  -ldl -le-hal
ifeq ($(INT_LIB),mpir)
       MPIR_DEF=-DUSE_MPIR
       CFLAGS +=$(MPIR_DEF)
       CXXFLAGS +=$(MPIR_DEF)
       LIBS+=-lmpir
else
       LIBS+=-lgmp -lgmpxx
endif
ifeq ($(OSVERSION),Linux)
	LIBS += -lrt
	ARCHFLAGS = -marm -march=armv7-a -mfpu=neon -mtune=cortex-a9
#	ARCHFLAGS = -march=corei7-avx
#	ARCHFLAGS = -march=corei7 -mno-avx
endif

BUILD_ARCH ?= $(ARCHFLAGS)

CFLAGS += $(BUILD_ARCH)
CXXFLAGS += $(BUILD_ARCH)

ifeq ($(OSVERSION),FreeBSD)
	CXX = clang++
	CC = clang
	CFLAGS += -DHAVE_DECL_LE32DEC -march=native
	CXXFLAGS += -DHAVE_DECL_LE32DEC -march=native
endif

# You might need to edit these paths too
LIBPATHS = -L/usr/local/lib -L/usr/lib -L/DBA/openssl/1.0.1f/lib/ -L/opt/adapteva/esdk/tools/host/lib
INCLUDEPATHS = -I/usr/local/include -I/usr/include -IxptMiner/includes/ -Iepiphany/src -I/opt/adapteva/esdk/tools/host/include

ifeq ($(OSVERSION),Darwin)
	EXTENSION = -mac
	GOT_MACPORTS := $(shell which port)
ifdef GOT_MACPORTS
	LIBPATHS += -L/opt/local/lib
	INCLUDEPATHS += -I/opt/local/include
endif
else
       EXTENSION =

endif

JHLIB = xptMiner/jhlib.o \

OBJS = \
        xptMiner/ticker.o \
	xptMiner/main.o \
	xptMiner/sha2.o \
	xptMiner/xptClient.o \
	xptMiner/xptClientPacketHandler.o \
	xptMiner/xptPacketbuffer.o \
	xptMiner/xptServer.o \
	xptMiner/xptServerPacketHandler.o \
	xptMiner/transaction.o \
	xptMiner/rh_riecoin.o \
	xptMiner/riecoinMiner.o


ifeq ($(ENABLE_OPENCL),1)
	OBJS += xptMiner/OpenCLObjects.o 
        OBJS += xptMiner/openCL.o
	CXXFLAGS += -DUSE_OPENCL
	LIBS += -lOpenCL
else

endif

all: xptminer$(EXTENSION) xptMiner/test epiphany/bin/e_primetest.elf

xptMiner/%.o: xptMiner/%.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDEPATHS) $< -o $@ 

xptMiner/rh_riecoin.o: xptMiner/rh_riecoin.c epiphany/src/common.h epiphany/src/ptest_data.h epiphany/src/modp_data.h
	$(CC) -c $(CFLAGS) $(INCLUDEPATHS) xptMiner/rh_riecoin.c -o $@ 

xptminer$(EXTENSION): $(OBJS:xptMiner/%=xptMiner/%) $(JHLIB:xptMiner/jhlib/%=xptMiner/jhlib/%)
	$(CXX) $(CFLAGS) $(LIBPATHS) $(INCLUDEPATHS) $(STATIC) -o $@ $^ $(LIBS) -flto

xptMiner/test: xptMiner/testharness.cpp xptMiner/riecoinMiner.o xptMiner/rh_riecoin.o xptMiner/sha2.o
	cd xptMiner && ./buildtest.sh

epiphany/bin/e_primetest.elf: epiphany/src/e_primetest.c epiphany/src/e_modp.c epiphany/src/common.h epiphany/src/ptest_data.h epiphany/src/modp_data.h epiphany/src/e_common.c
	cd epiphany && ./build.sh

clean:
	-rm -f xptminer
	-rm -f xptMiner/*.o
	-rm -f xptMiner/jhlib/*.o
