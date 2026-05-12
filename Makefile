CXX = g++
CATIMA_PREFIX ?= /usr/local
CFLAGS = $(shell root-config --cflags) -Isrc/ -I$(CATIMA_PREFIX)/include
LIBS = $(shell root-config --glibs) -lGeom -lEve -lRGL -L$(CATIMA_PREFIX)/lib -lcatima
OBJS = lib/main.o lib/MUSIC_Simulator.o lib/NuclideFinder.o lib/Particle.o lib/EnergyLoss.o lib/FourVector.o

all: musicsim

musicsim: $(OBJS) | lib
	$(CXX) -o $@ $(OBJS) $(LIBS)

lib:
	mkdir -p lib

lib/main.o: src/main.cpp src/MUSIC_Simulator.hpp | lib
	$(CXX) $(CFLAGS) -c src/main.cpp -o $@

lib/MUSIC_Simulator.o: src/MUSIC_Simulator.cpp src/MUSIC_Simulator.hpp src/NuclideFinder.hpp src/Particle.hpp | lib
	$(CXX) $(CFLAGS) -c src/MUSIC_Simulator.cpp -o $@

lib/Particle.o: src/Particle.cpp src/Particle.hpp src/FourVector.hpp src/EnergyLoss.hpp | lib
	$(CXX) $(CFLAGS) -c src/Particle.cpp -o $@

lib/FourVector.o: src/FourVector.cpp src/FourVector.hpp | lib
	$(CXX) -Isrc/ -c src/FourVector.cpp -o $@

lib/EnergyLoss.o: src/EnergyLoss.cpp src/EnergyLoss.hpp | lib
	$(CXX) $(CFLAGS) -c src/EnergyLoss.cpp -o $@

lib/NuclideFinder.o: src/NuclideFinder.cpp src/NuclideFinder.hpp | lib
	$(CXX) -Isrc/ -c src/NuclideFinder.cpp -o $@

clean:
	-rm -f musicsim lib/*.o
