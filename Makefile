# Makefile for the Argonne Nuclear Data Exploration Software (ANDES)
# Created by Daniel Santiago-Gonzalez (dsg@anl.gov) Nov/2020
# Use 'make' to generate the andes executable file.
# Use 'make clean' to rm andes and associated objects (*.o)
CXX = g++
CFLAGS = $(shell root-config --cflags) -Isrc/
LIBS = $(shell root-config --glibs) -lEve -lRGL
OBJS = lib/main.o lib/MUSIC_Simulator.o lib/NuclideFinder.o

musicsim: $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LIBS)

lib/main.o: src/main.cpp lib/MUSIC_Simulator.o lib/NuclideFinder.o
	$(CXX) $(CFLAGS) -c src/main.cpp -o $@

lib/MUSIC_Simulator.o: src/MUSIC_Simulator.cpp lib/Particle.o lib/NuclideFinder.o
	$(CXX) $(CFLAGS) -c src/MUSIC_Simulator.cpp -o $@

lib/Particle.o: src/Particle.cpp lib/FourVector.o lib/EnergyLoss.o
	$(CXX) $(CFLAGS) -c src/Particle.cpp -o $@

lib/FourVector.o: src/FourVector.cpp src/FourVector.hpp
	$(CXX) -Isrc/ -c src/FourVector.cpp -o $@

lib/EnergyLoss.o: src/EnergyLoss.cpp src/EnergyLoss.hpp
	$(CXX) $(CFLAGS) -c src/EnergyLoss.cpp -o $@

lib/NuclideFinder.o: src/NuclideFinder.cpp src/NuclideFinder.hpp
	$(CXX) -Isrc/ -c src/NuclideFinder.cpp -o $@


clean:
	-rm andes lib/*.o
