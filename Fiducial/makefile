#ROOT_PATH = $ROOTSYS
ROOT_CFLAGS =$(shell root-config --cflags)
ROOT_LIBS = $(shell root-config --glibs)

FIDUCIAL_INTERFACE = interface


all:

	g++ -g $(ROOT_CFLAGS) -I $(FIDUCIAL_INTERFACE) $(ROOT_LIBS)  -o test/CommonFiducialSkim.exe src/common_fiducial_skim.cc src/object_cuts.cc src/event_cuts.cc src/histogram_builder.cc

clean:
	rm src/CommonFiducialSkim.exe

#g++ -g -Wall -fPIC -pthread -m64  $(ROOT_CFLAGS) -I $(FIDUCIAL_INTERFACE) $(ROOT_LIBS)  -o test/CommonFiducialSkim.exe src/common_fiducial_skim.cc src/object_cuts.cc src/event_cuts.cc src/histogram_builder.cc

#g++ src/common_fiducial_skim.cc src/mc_particle_identification.cc $(ROOT_CFLAGS) $(ROOT_LIBS) -I $(FIDUCIAL_INTERFACE) -o test/CommonFiducialSkim.exe

#g++ src/test.cpp $(ROOT_CFLAGS) $(ROOT_LIBS) -I $(ROOT_PATH) -o test/test.exe