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

# srim-cache is a separate tool with its own main(); keep it out of the
# simulator objects or the link sees two.
SRIMCACHE_SRC = $(SRCDIR)/srim-cache.cpp
SOURCES = $(filter-out $(SRIMCACHE_SRC),$(wildcard $(SRCDIR)/*.cpp))
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
HEADERS = $(wildcard $(INCDIR)/*.hpp)

all: musicsim srim-cache

musicsim: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LIBS)

srim-cache: $(SRIMCACHE_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS) | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) musicsim srim-cache

.PHONY: all clean
