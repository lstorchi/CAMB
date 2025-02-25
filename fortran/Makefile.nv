#CAMB Makefile

#Set FISHER=Y to compile bispectrum fisher matrix code
FISHER=

#Set FORUTILSPATH to the path where the libforutils.a file can be found.
#The OUTPUT_DIR will be appended.

ifneq "$(wildcard ../forutils)" ""
FORUTILSPATH ?= $(shell pwd)/../forutils
else
ifneq "$(wildcard ../../forutils)" ""
FORUTILSPATH ?= $(shell pwd)/../../forutils
else
ifneq "$(wildcard ../../../forutils)" ""
FORUTILSPATH ?= $(shell pwd)/../../../forutils
endif
endif
endif

ifeq ($(FORUTILSPATH),)
$(error Use  git clone --recurse-submodules, install forutils from https://github.com/cmbant/forutils, or set FORUTILSPATH variable)
endif

#native optimization does not work on Mac gfortran or heterogeneous clusters
CLUSTER_SAFE ?= 0
ifneq ($(CLUSTER_SAFE), 0)
NONNATIVE = 1
endif

NONNATIVE = 1

COMMON_FFLAGS = -cpp -Minform=warn -Mextend
FFLAGS = -O0 -g $(COMMON_FFLAGS)
DEBUGFLAGS = -g -fbacktrace -ffpe-trap=invalid,overflow,zero -fbounds-check $(COMMON_FFLAGS)
MODOUT =  -I$(OUTPUT_DIR)
SMODOUT = -J$(DLL_DIR)

ifneq ($(FISHER),)
F90CRLINK += -lblas -llapack
endif

ifndef NONNATIVE
#Note this seems to make code slightly slower in some cases, use CLUSTER_SAFE=1 to test without
FFLAGS+=-march=native
endif

IFLAG = -I

#Settings for building camb_fits
#Location of FITSIO and name of library
FITSDIR       ?= /usr/local/lib
FITSLIB       = cfitsio
#Location of HEALPIX for building camb_fits
HEALPIXDIR    ?= /usr/local/healpix

ifneq ($(FISHER),)
# Its dependencies are all meet by the libutils.a which always added.
FFLAGS += -DFISHER
endif

DEBUGFLAGS ?= $(FFLAGS)
Debug: FFLAGS = $(DEBUGFLAGS)

COMPILER = nvfortran
F90C     = nvfortran

include ./Makefile_main
