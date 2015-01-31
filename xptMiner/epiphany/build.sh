#!/bin/bash

set -e

ESDK=${EPIPHANY_HOME}
ELIBS=${ESDK}/tools/host/lib
EINCS=${ESDK}/tools/host/include
ELDF=${ESDK}/bsps/current/fast.ldf

P_TO_TEST=49979693

# Create the binaries directory
mkdir -p bin/

CROSS_PREFIX=
case $(uname -p) in
	arm*)
		# Use native arm compiler (no cross prefix)
		CROSS_PREFIX=
		;;
	   *)
		# Use cross compiler
		CROSS_PREFIX="arm-linux-gnueabihf-"
		;;
esac

# Build HOST side application
#${CROSS_PREFIX}gcc -O3 -std=gnu99 src/modptest.c -o bin/modptest.elf -I ${EINCS} -L ${ELIBS} -le-hal  -lm -lpthread -DP_TO_TEST=${P_TO_TEST}

# Build DEVICE side program
e-gcc -O3 -g -T ${ELDF} -std=gnu99 src/e_modp.c -o bin/e_modp.elf -le-lib -lm -ffast-math -Wall -mfp-mode=int -DNO_RESULT_DEBUG
#e-gcc -O3 -S -T ${ELDF} -std=gnu99 src/e_modp.c -o bin/e_modp.s -le-lib -lm -ffast-math -mfp-mode=int
e-gcc -O3 -g -T ${ELDF} -std=gnu99 src/e_primetest.c -o bin/e_primetest.elf -le-lib -lm -ffast-math -Wall -mfp-mode=int -DNO_RESULT_DEBUG
e-gcc -O3 -S -T ${ELDF} -std=gnu99 src/e_primetest.c -o bin/e_primetest.s -le-lib -lm -ffast-math -mfp-mode=int

# Convert ebinary to SREC file
e-objcopy --srec-forceS3 --output-target srec bin/e_modp.elf bin/e_modp.srec
e-objcopy --srec-forceS3 --output-target srec bin/e_primetest.elf bin/e_primetest.srec


