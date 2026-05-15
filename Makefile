CXX           = g++
CATIMA_PREFIX ?= /usr/local
VERSION       ?= dev

SRCDIR = src
INCDIR = include
OBJDIR = lib

CXXFLAGS = $(shell root-config --cflags) -I$(INCDIR) -I$(CATIMA_PREFIX)/include \
           -DMUSICSIM_VERSION=\"$(VERSION)\"
LIBS     = $(shell root-config --glibs) -lGeom -lEve -lRGL -lMathMore \
           -L$(CATIMA_PREFIX)/lib -lcatima

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
HEADERS = $(wildcard $(INCDIR)/*.hpp)

all: musicsim

musicsim: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS) | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) musicsim

.PHONY: all clean
