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

  setenv DEVELOPER marc
  setenv GNU_MAKE make
  setenv TARGET shared,opt

  setenv DV_DIR ${CCSE_DIR}/Combustion/Chemistry/tools/fuego/Pythia
  setenv BLD_ROOT ${DV_DIR}/builds
  setenv BLD_CONFIG ${DV_DIR}/config
  setenv PYTHIA_DIR ${DV_DIR}/pythia
  export TMP_FUEGO_DIR=`pwd`/Fuego

  setenv EXPORT_ROOT ${DV_DIR}/pythia-0.4

  set path=( ${BLD_CONFIG}/make ${path} ${EXPORT_ROOT}/bin )
#  setenv PYTHONPATH ${EXPORT_ROOT}/modules
  setenv PYTHONPATH ${TMP_FUEGO_DIR}:${EXPORT_ROOT}/packages/fuego:${EXPORT_ROOT}/packages/journal:${EXPORT_ROOT}/packages/pyre:${EXPORT_ROOT}/packages/weaver

  setenv LD_LIBRARY_PATH ${EXPORT_ROOT}/lib:${LD_LIBRARY_PATH}

# Python

  setenv PYTHON_DIR /usr
  setenv PYTHON_VERSION 2.5
  setenv PYTHON_INCDIR ${PYTHON_DIR}/include/python${PYTHON_VERSION}
  setenv PYTHON_LIBDIR ${PYTHON_DIR}/lib/python${PYTHON_VERSION}

  rehash

# version
# $Id$

# Generated automatically by ShMill on Sat May 17 15:06:54 2003

# End of file 
