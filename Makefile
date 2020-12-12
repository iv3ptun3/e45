# $Id: GNUmakefile,v 1.1 1999/01/07 16:05:40 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name		:= hyptpc1
G4TARGET	:= $(name)
G4EXLIB		:= true
# CPPVERBOSE	:= true

PHONY: all
all: lib bin

#ifndef G4WORKDIR
#  G4WORKDIR	= ./
#endif
  G4WORKDIR	= ./

include $(G4INSTALL)/config/binmake.gmk

SOFLAGS		+= -shared
ROOTLIBS	:= $(shell root-config --libs)
ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTGLIBS	:= $(shell root-config --glibs)
# CXXFLAGS	+= -DDEBUG
CXXFLAGS	+= $(ROOTCFLAGS)
CPPFLAGS	+= $(ROOTCFLAGS)

# CERN_ROOT = /sw
# CLIB = -L$(CERN_ROOT)/lib -lpawlib -lgraflib -lgrafX11 \
#         -lpacklib -lphtools -lmathlib -lkernlib -lnsl
#XLIB = -L/usr/lib -lXt  -lX11 -lXp -lXext -lm -lc

FC = gfortran -m32
FFLAGS +=  -c ./include

LDLIBS   += $(ROOTLIBS)

##LDLIBS +=$(CLIB)
#ifndef G4INSTALL
#  G4INSTALL = ../
#endif
