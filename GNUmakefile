export G4LIB_BUILD_STATIC=1
export G4LIB_BUILD_SHARED=

CPPFLAGS += -I$(ROOTSYS)/include 
EXTRALIBS = $(shell root-config --glibs)

name := MINER-sim
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

CPPFLAGS += -no-pie -I$(CRYHOME)/src
EXTRALIBS  += -fPIC -L$(CRYHOME)/lib -lCRY

# CPPFLAGS += -I$(ROOTSYS)/include -Wl,--no-as-needed 
# EXTRALIBS = $(shell root-config --glibs)

.PHONY: all
all: lib bin

include $(G4INSTALL)/config/binmake.gmk
