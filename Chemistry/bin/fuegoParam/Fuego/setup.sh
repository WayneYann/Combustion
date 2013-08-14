#!/bin/sh
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Michael A.G. Aivazis
#                      California Institute of Technology
#                      (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

# build system

  export DEVELOPER=marc
  export GNU_MAKE=make
  export TARGET=shared,opt

  export DV_DIR=${CCSE_DIR}/Combustion/Chemistry/tools/fuego/Pythia
  export BLD_ROOT=${DV_DIR}/builds
  export BLD_CONFIG=${DV_DIR}/config
  export PYTHIA_DIR=${DV_DIR}/pythia
  export TMP_FUEGO_DIR=`pwd`/Fuego

  export EXPORT_ROOT=${DV_DIR}/pythia-0.4

  export PATH=$BLD_CONFIG/make:${PATH}:${EXPORT_ROOT}/bin
  export PYTHONPATH=${TMP_FUEGO_DIR}:${EXPORT_ROOT}/packages/fuego:${EXPORT_ROOT}/packages/journal:${EXPORT_ROOT}/packages/pyre:${EXPORT_ROOT}/packages/weaver
  export LD_LIBRARY_PATH=${EXPORT_ROOT}/lib:${LD_LIBRARY_PATH}

# Python

  export PYTHON_DIR=/usr
  export PYTHON_VERSION=2.6
  export PYTHON_INCDIR=${PYTHON_DIR}/include/python${PYTHON_VERSION}
  export PYTHON_LIBDIR=${PYTHON_DIR}/lib/python${PYTHON_VERSION}


# version
# $Id$

# Generated automatically by ShMill on Sat May 17 15:06:54 2003

# End of file 
