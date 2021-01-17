# Makefile for the MUSIC Simulator
# Created by Daniel Santiago-Gonzalez (dsg@anl.gov) Jan/2021
# Use 'make musicsim' to generate the musicsim executable file.
# Use 'make make-srim-table' to generate the make-srim-table executable file.
# Use 'make' or 'make all' to generate both executable files (default)
# Use 'make clean' to rm musicsim, makes-srim-table, and associated objects (*.o)
CXX = g++
CFLAGS = $(shell root-config --cflags) -Isrc/
LIBS = $(shell root-config --glibs) -lGeom -lEve -lRGL 
OBJS = lib/main.o lib/MUSIC_Simulator.o lib/NuclideFinder.o lib/Particle.o lib/EnergyLoss.o lib/FourVector.o
#OBJS = lib/main.o lib/MUSIC_Simulator.o 

all: musicsim make-srim-table

musicsim: $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LIBS)

make-srim-table: lib/make-srim-table.o lib/SRIM_Table_Maker.o lib/NuclideFinder.o
	$(CXX) -o $@ lib/make-srim-table.o lib/SRIM_Table_Maker.o lib/NuclideFinder.o

lib/make-srim-table.o: src/make-srim-table.cpp lib/SRIM_Table_Maker.o lib/NuclideFinder.o 
	$(CXX) $(CFLAGS) -c src/make-srim-table.cpp -o $@

lib/SRIM_Table_Maker.o: src/SRIM_Table_Maker.cpp lib/NuclideFinder.o
	$(CXX) -Isrc/ -c src/SRIM_Table_Maker.cpp -o $@

lib/main.o: src/main.cpp lib/MUSIC_Simulator.o lib/NuclideFinder.o 
	$(CXX) $(CFLAGS) -c src/main.cpp -o $@

lib/MUSIC_Simulator.o: src/MUSIC_Simulator.cpp lib/NuclideFinder.o lib/Particle.o 
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
	-rm musicsim make-srim-table lib/*.o
