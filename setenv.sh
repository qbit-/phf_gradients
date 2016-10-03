#!/bin/sh
#. /opt/pgi/linux86-64/13.10/pgi.sh
gdvroot=/home/shiva/software/gdv-h21
GAUSS_MEMDEF=67108864
GAUSS_SCRDIR=/tmp
export PATH /opt/pgi/linux86-64/13.3/bin:$PATH
export gdvroot GAUSS_MEMDEF GAUSS_SCRDIR
. $gdvroot/gdv/bsd/gdv.profile
