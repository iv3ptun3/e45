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

SOFLAGS		:= -shared

# changed by hjb
ROOTCFLAGS := $(shell root-config --cflags)
ROOTINCLUDE := -I/opt/homebrew/Cellar/root/6.32.08/include/root
ROOTGLIBS    := $(shell root-config --glibs)

CXX         := clang++
CXXFLAGS    += -std=c++17 -stdlib=libc++ -pthread -arch arm64
CXXFLAGS    += -Wall -Wextra
CPPFLAGS    += $(ROOTINCLUDE)

LDFLAGS     += -L/opt/homebrew/Cellar/root/6.32.08/lib/root \
           -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint \
           -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMinuit
# LDFLAGS     += -L/opt/homebrew/Cellar/root/6.32.08/lib/root -lMinuit2
LDFLAGS     += $(ROOTGLIBS)

#ROOTLIBS	:= $(shell root-config --libs)

#ROOTCFLAGS	:= $(shell root-config --cflags)
#ROOTGLIBS	:= $(shell root-config --glibs)
# CXXFLAGS	+= -DDEBUG
#CXXFLAGS	+= $(ROOTCFLAGS)
#CPPFLAGS	+= $(ROOTCFLAGS)

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
