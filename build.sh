#!/bin/bash

#  Modified from build.sh by:
#  Jacob Williams : 2/8/2014
#     - modified 6/23/2014

SRCDIR='src/'                #source directory
CMSRC='/Users/doughunsaker/Documents/XCodeProjects/CommonSource'       #common source directory
BUILDDIR='obj/'              #build directory for library
BINDIR='bin/'                #build directory for executable
FEXT='.f90'                  #fortran file extension
OBJEXT='.o'                  #object code extension
LIBEXT='.a'                  #static library extension
MODEXT='.mod'                #fortran module file extension
WC='*'                       #wildcard character
EXEOUT='MachUp.out'              #name of output program

#
# Always a clean build:
#

mkdir -p $BINDIR
mkdir -p $BUILDDIR

rm -f $BUILDDIR$WC$OBJEXT
rm -f $BUILDDIR$WC$MODEXT
rm -f $BUILDDIR$WC$LIBEXT

cd $SRCDIR
gfortran -fdefault-real-8 $CMSRC/json.f90 $CMSRC/myjson.f90 $CMSRC/math.f90 $CMSRC/dataset.f90 $CMSRC/airfoil.f90 $CMSRC/section.f90 wing.f90 plane.f90 view.f90 special_functions.f90 loads.f90 main.f90 -o $EXEOUT
cd ..
mv $SRCDIR$WC$MODEXT $BUILDDIR
mv $SRCDIR$EXEOUT $BINDIR

