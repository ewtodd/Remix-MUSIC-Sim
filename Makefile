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

# Command-line tools have their own main(); keep them out of the simulator
# objects or the link sees duplicate entry points.
SRIMCACHE_SRC = $(SRCDIR)/srim-cache.cpp
CONVERTER_SRC = $(SRCDIR)/legacy-msc-to-toml.cpp
TOOL_SOURCES = $(SRIMCACHE_SRC) $(CONVERTER_SRC)
SOURCES = $(filter-out $(TOOL_SOURCES),$(wildcard $(SRCDIR)/*.cpp))
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
HEADERS = $(wildcard $(INCDIR)/*.hpp)
TESTDIR = $(OBJDIR)/tests
ROOT_CHECK = $(TESTDIR)/root_checks
VAVILOV_CHECK = $(TESTDIR)/vavilov_statistics
NUCLIDE_CHECK = $(TESTDIR)/nuclide_loader
SANITIZE_DIR = $(TESTDIR)/sanitize
SANITIZE_OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(SANITIZE_DIR)/%.o)
# GCC 15's sanitizer passes are pathologically slow on the generated nuclide
# table at -O1. Instrumented -O0 builds retain the checks and finish quickly;
# Nix injects _FORTIFY_SOURCE, whose expected -O0 preprocessor warning is muted.
SANITIZE_FLAGS = -O0 -g1 -Wno-cpp -fsanitize=address,undefined \
                 -fno-omit-frame-pointer
SANITIZE_SIM = $(SANITIZE_DIR)/musicsim
SANITIZE_CONVERTER = $(SANITIZE_DIR)/legacy-msc-to-toml
SANITIZE_ROOT_CHECK = $(SANITIZE_DIR)/root_checks
SANITIZE_VAVILOV_CHECK = $(SANITIZE_DIR)/vavilov_statistics
SANITIZE_NUCLIDE_CHECK = $(SANITIZE_DIR)/nuclide_loader
CPP_FORMAT_FILES = $(sort $(wildcard $(INCDIR)/*.hpp $(SRCDIR)/*.cpp tests/*.cpp))
TOML_FORMAT_FILES = $(sort basic.toml $(wildcard ControlExamples/*/*.toml tests/controls/*.toml))
CLANG_FORMAT_STYLE = {BasedOnStyle: LLVM, BreakStringLiterals: false}
STRICT_FLAGS = -Wall -Wextra -Wpedantic -Wconversion -Wshadow \
               -Wnon-virtual-dtor -Wold-style-cast -Woverloaded-virtual \
               -Wnull-dereference -Wdouble-promotion -Wformat=2 -Werror

all: musicsim srim-cache legacy-msc-to-toml

musicsim: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LIBS)

# The SRIM table generator srim-cache drives: nix passes SRIM-nix's
# make-srim-table store path; empty means whatever is on PATH at run time.
SRIM_TABLE_BIN ?=

srim-cache: $(SRIMCACHE_SRC)
	$(CXX) $(CXXFLAGS) -DMUSICSIM_SRIM_TABLE_BIN='"$(SRIM_TABLE_BIN)"' -o $@ $<

legacy-msc-to-toml: $(CONVERTER_SRC)
	$(CXX) -std=c++17 -Wall -Wextra -Wpedantic -o $@ $<

$(ROOT_CHECK): tests/root_checks.cpp | $(TESTDIR)
	$(CXX) -std=c++17 -Wall -Wextra -Wpedantic $(shell root-config --cflags) \
		-o $@ $< $(shell root-config --libs)

$(VAVILOV_CHECK): tests/vavilov_statistics.cpp src/VavilovSampler.cpp \
		include/VavilovSampler.hpp | $(TESTDIR)
	$(CXX) -std=c++17 -Wall -Wextra -Wpedantic $(shell root-config --cflags) \
		-I$(INCDIR) -o $@ tests/vavilov_statistics.cpp src/VavilovSampler.cpp \
		$(shell root-config --libs) -lMathMore

$(NUCLIDE_CHECK): tests/nuclide_loader.cpp src/NuclideFinder.cpp \
		include/NuclideFinder.hpp | $(TESTDIR)
	$(CXX) -std=c++17 -Wall -Wextra -Wpedantic $(shell root-config --cflags) \
		-I$(INCDIR) -o $@ tests/nuclide_loader.cpp src/NuclideFinder.cpp \
		$(shell root-config --libs)

$(TESTDIR):
	mkdir -p $(TESTDIR)

# Build-only target for the check binaries; the nix build invokes this so
# `nix build` ships them next to musicsim and the suite can run without make.
test-bins: $(ROOT_CHECK) $(VAVILOV_CHECK) $(NUCLIDE_CHECK)

check: all $(ROOT_CHECK) $(VAVILOV_CHECK) $(NUCLIDE_CHECK)
	MUSICSIM_BIN=$(abspath musicsim) \
	CONVERTER_BIN=$(abspath legacy-msc-to-toml) \
	ROOT_CHECK_BIN=$(abspath $(ROOT_CHECK)) \
	VAVILOV_CHECK_BIN=$(abspath $(VAVILOV_CHECK)) \
	NUCLIDE_CHECK_BIN=$(abspath $(NUCLIDE_CHECK)) \
		bash tests/run-tests.sh

strict:
	@command -v $(CXX) >/dev/null
	@for source in $(SOURCES) $(TOOL_SOURCES) tests/*.cpp; do \
		echo "strict: $$source"; \
		$(CXX) $(CXXFLAGS) $(STRICT_FLAGS) -fsyntax-only $$source || exit 1; \
	done

$(SANITIZE_DIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS) | $(SANITIZE_DIR)
	$(CXX) $(CXXFLAGS) $(STRICT_FLAGS) $(SANITIZE_FLAGS) -c $< -o $@

$(SANITIZE_SIM): $(SANITIZE_OBJECTS)
	$(CXX) $(SANITIZE_FLAGS) -o $@ $^ $(LIBS)

$(SANITIZE_CONVERTER): $(CONVERTER_SRC) | $(SANITIZE_DIR)
	$(CXX) -std=c++17 $(STRICT_FLAGS) $(SANITIZE_FLAGS) -o $@ $<

$(SANITIZE_ROOT_CHECK): tests/root_checks.cpp | $(SANITIZE_DIR)
	$(CXX) -std=c++17 $(STRICT_FLAGS) $(SANITIZE_FLAGS) \
		$(shell root-config --cflags) -o $@ $< $(shell root-config --libs)

$(SANITIZE_VAVILOV_CHECK): tests/vavilov_statistics.cpp \
		$(SANITIZE_DIR)/VavilovSampler.o include/VavilovSampler.hpp | $(SANITIZE_DIR)
	$(CXX) -std=c++17 $(STRICT_FLAGS) $(SANITIZE_FLAGS) \
		$(shell root-config --cflags) -I$(INCDIR) -o $@ \
		tests/vavilov_statistics.cpp $(SANITIZE_DIR)/VavilovSampler.o \
		$(shell root-config --libs) -lMathMore

$(SANITIZE_NUCLIDE_CHECK): tests/nuclide_loader.cpp \
		$(SANITIZE_DIR)/NuclideFinder.o \
		include/NuclideFinder.hpp | $(SANITIZE_DIR)
	$(CXX) -std=c++17 $(STRICT_FLAGS) $(SANITIZE_FLAGS) \
		$(shell root-config --cflags) -I$(INCDIR) -o $@ \
		tests/nuclide_loader.cpp $(SANITIZE_DIR)/NuclideFinder.o \
		$(shell root-config --libs)

$(SANITIZE_DIR):
	mkdir -p $(SANITIZE_DIR)

sanitize-check: $(SANITIZE_SIM) $(SANITIZE_CONVERTER) \
		$(SANITIZE_ROOT_CHECK) $(SANITIZE_VAVILOV_CHECK) \
		$(SANITIZE_NUCLIDE_CHECK)
	ASAN_OPTIONS=detect_leaks=0:halt_on_error=1 \
	UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
	MUSICSIM_BIN=$(abspath $(SANITIZE_SIM)) \
	CONVERTER_BIN=$(abspath $(SANITIZE_CONVERTER)) \
	ROOT_CHECK_BIN=$(abspath $(SANITIZE_ROOT_CHECK)) \
	VAVILOV_CHECK_BIN=$(abspath $(SANITIZE_VAVILOV_CHECK)) \
	NUCLIDE_CHECK_BIN=$(abspath $(SANITIZE_NUCLIDE_CHECK)) \
		bash tests/run-tests.sh

format:
	@command -v clang-format >/dev/null || { \
		echo "clang-format is required (enter 'nix develop')" >&2; exit 1; }
	@command -v taplo >/dev/null || { \
		echo "taplo is required (enter 'nix develop')" >&2; exit 1; }
	clang-format -i --style='$(CLANG_FORMAT_STYLE)' $(CPP_FORMAT_FILES)
	taplo fmt $(TOML_FORMAT_FILES)

format-check:
	@command -v clang-format >/dev/null || { \
		echo "clang-format is required (enter 'nix develop')" >&2; exit 1; }
	@command -v taplo >/dev/null || { \
		echo "taplo is required (enter 'nix develop')" >&2; exit 1; }
	clang-format --dry-run --Werror --style='$(CLANG_FORMAT_STYLE)' \
		$(CPP_FORMAT_FILES)
	taplo check $(TOML_FORMAT_FILES)
	taplo fmt --check $(TOML_FORMAT_FILES)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS) | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf $(OBJDIR) musicsim srim-cache legacy-msc-to-toml

.PHONY: all check test-bins strict sanitize-check format format-check clean
