#!/bin/bash

set -e

ESDK=${EPIPHANY_HOME}
if [ .$ESDK == . ]; then ESDK=/opt/adapteva/esdk; fi
ELIBS=${ESDK}/tools/host/lib:${LD_LIBRARY_PATH}
EHDF=${EPIPHANY_HDF}
if [ .$EHDF == . ]; then EHDF=/opt/adapteva/esdk/bsps/current/platform.hdf; fi

sudo -E LD_LIBRARY_PATH=${ELIBS} EPIPHANY_HDF=${EHDF} stdbuf -o L ./xptminer -o http://mining.ypool.net -t 1 $* | tee -a run.log

