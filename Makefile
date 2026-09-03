CXX           = g++
CATIMA_PREFIX ?= /usr/local
VERSION       ?= dev

SRCDIR = src
INCDIR = include
OBJDIR = lib

# Where the wordmark (assets/remix.txt) lives. Baked in at build so the
# binary stays self-locating regardless of the run-time cwd; nix points it
# at the installed copy.
ASSETS_DIR     ?= $(abspath assets)
ASSETS_DIR_OUT ?= $(ASSETS_DIR)

CXXFLAGS = $(shell root-config --cflags) -I$(INCDIR) -I$(CATIMA_PREFIX)/include \
           -DMUSICSIM_VERSION=\"$(VERSION)\" \
           -DMUSICSIM_ASSETS_DIR='"$(ASSETS_DIR_OUT)"'
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

# The SRIM table generator srim-cache drives: nix passes SRIM-nix's
# make-srim-table store path; empty means whatever is on PATH at run time.
SRIM_TABLE_BIN ?=

srim-cache: $(SRIMCACHE_SRC)
	$(CXX) $(CXXFLAGS) -DMUSICSIM_SRIM_TABLE_BIN='"$(SRIM_TABLE_BIN)"' -o $@ $<

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS) | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) musicsim srim-cache

.PHONY: all clean
