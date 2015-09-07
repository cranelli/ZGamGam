#ROOT_PATH = $ROOTSYS
ROOT_CFLAGS =$(shell root-config --cflags)
ROOT_LIBS = $(shell root-config --glibs)

FIDUCIAL_INTERFACE = interface


all: common dress gen_dress

#	g++ -g $(ROOT_CFLAGS) -I $(FIDUCIAL_INTERFACE) $(ROOT_LIBS)  -o test/CommonFiducialSkim.exe src/common_fiducial_skim.cc src/object_cuts.cc src/event_cuts.cc src/histogram_builder.cc

common:
	g++ -g $(ROOT_CFLAGS) -I $(FIDUCIAL_INTERFACE) $(ROOT_LIBS)  -o test/CommonFiducialSkim.exe src/common_fiducial_skim.cc src/object_cuts.cc src/event_cuts.cc src/histogram_builder.cc

dress:
	g++ -g $(ROOT_CFLAGS) -I $(FIDUCIAL_INTERFACE) $(ROOT_LIBS)  -o test/CommonFiducialDressLeptonSkim.exe src/common_fiducial_dress_lepton_skim.cc src/object_cuts.cc src/event_cuts.cc src/histogram_builder.cc

gen:
	g++ -g $(ROOT_CFLAGS) -I $(FIDUCIAL_INTERFACE) $(ROOT_LIBS)  -o test/GenHistograms.exe src/gen_histograms.cc src/object_cuts.cc src/event_cuts.cc src/histogram_builder.cc

gen_dress:
	g++ -g $(ROOT_CFLAGS) -I $(FIDUCIAL_INTERFACE) $(ROOT_LIBS)  -o test/MakeGenDressHistograms.exe src/make_gen_dress_histograms.cc src/object_cuts.cc src/event_cuts.cc src/histogram_builder.cc

clean:
	rm test/CommonFiducialSkim.exe
	rm test/CommonFiducialDressLeptonSkim.exe

#g++ -g -Wall -fPIC -pthread -m64  $(ROOT_CFLAGS) -I $(FIDUCIAL_INTERFACE) $(ROOT_LIBS)  -o test/CommonFiducialSkim.exe src/common_fiducial_skim.cc src/object_cuts.cc src/event_cuts.cc src/histogram_builder.cc

#g++ src/common_fiducial_skim.cc src/mc_particle_identification.cc $(ROOT_CFLAGS) $(ROOT_LIBS) -I $(FIDUCIAL_INTERFACE) -o test/CommonFiducialSkim.exe

#g++ src/test.cpp $(ROOT_CFLAGS) $(ROOT_LIBS) -I $(ROOT_PATH) -o test/test.exe