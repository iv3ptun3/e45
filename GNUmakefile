# $Id: GNUmakefile,v 1.1 1999/01/07 16:05:40 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := hyptpc1
G4TARGET := $(name)
G4EXLIB := true

ifndef G4WORKDIR
  G4WORKDIR = ./
endif

PHONY: all
all: lib bin


G4WORKDIR = ./
SOFLAGS        += -shared
ROOTLIBS       := $(shell root-config --libs)
ROOTCFLAGS     := $(shell root-config --cflags)
ROOTGLIBS      := $(shell root-config --glibs)
CXXFLAGS        = -O -Wall -fPIC
CXXFLAGS       += $(ROOTCFLAGS)
CPPFLAGS       += $(ROOTCFLAGS)

#include $(G4INSTALL)/config/binmake.gmk

ifndef G4INSTALL
#  G4INSTALL = /sw/packages/Geant4/geant4.9.4.p03
#  G4INSTALL = /sw/packages/geant4/9.6.4/share/Geant4-9.6.4/geant4make
  G4INSTALL = /sw/packages/geant4/9.5.2/share/Geant4-9.5.2/geant4make	
#  G4INSTALL = /sw/packages/geant4/10.2.2/share/Geant4-10.2.2/geant4make
endif


include $(G4INSTALL)/config/binmake.gmk

CERN_ROOT = /sw
CLIB = -L$(CERN_ROOT)/lib -lpawlib -lgraflib -lgrafX11 \
        -lpacklib -lphtools -lmathlib -lkernlib -lnsl
XLIB = -L/usr/lib -lXt  -lX11 -lXp -lXext -lm -lc

FC = gfortran -m32
FLAGS = -O -fno-automatic -ffixed-line-length-132
FFLAGS +=  -c ./include

CXXFLAGS += $(ROOTCFLAGS)
CPPFLAGS += -I$(ROOTSYS)/include 
LDLIBS   += $(ROOTLIBS) 
CPPFLAGS += -I/sw/include
LDLIBS   += -L/sw/lib


##LDLIBS +=$(CLIB)
#ifndef G4INSTALL
#  G4INSTALL = ../
#endif



