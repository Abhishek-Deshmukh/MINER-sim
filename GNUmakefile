name := MinerSim
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../../..
endif

.PHONY: all
all: lib bin

CPPFLAGS += -I$(ROOTSYS)/include
EXTRALIBS = $(shell root-config --glibs)

include $(G4INSTALL)/config/binmake.gmk