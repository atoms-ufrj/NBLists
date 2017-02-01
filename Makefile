# Define DEBUG or FAST mode:
#   Build the "fast" version with: `make` or `make DEBUG=0`
#   Build the "debug" version with: `make DEBUG=1`
DEBUG ?= 0

# Installation prefix:
PREFIX ?= /usr/local

# Compilers and their basic options:
FORT ?= gfortran
CC ?= gcc

BASIC_F_OPTS = -march=native -m64 -fPIC -fopenmp -cpp -fmax-errors=1
BASIC_C_OPTS = -march=native -m64 -fPIC -fopenmp -cpp -fmax-errors=1

# Warning-related options:
BASIC_F_OPTS += -Wall -Wno-maybe-uninitialized
BASIC_C_OPTS += -Wall -Wno-maybe-uninitialized

# Option FAST (default):
FAST_F_OPTS = -Ofast
FAST_C_OPTS = -Ofast

# Option DEBUG:
DEBUG_F_OPTS = --coverage -g -Og -fcheck=all -Ddebug
DEBUG_C_OPTS = --coverage -g -Og -fstack-check -fbounds-check -Ddebug

# Checks chosen option:
ifeq ($(DEBUG), 1)
F_OPTS = $(BASIC_F_OPTS) $(DEBUG_F_OPTS)
C_OPTS = $(BASIC_C_OPTS) $(DEBUG_C_OPTS)
else
F_OPTS = $(BASIC_F_OPTS) $(FAST_F_OPTS)
C_OPTS = $(BASIC_C_OPTS) $(FAST_C_OPTS)
endif

LN_INC_OPT = -I$(INCDIR)
LN_SO_OPT = -Wl,-rpath,'$$ORIGIN/../lib'
LN_OPTS = $(LN_INC_OPT) $(LN_SO_OPT)

SRCDIR = ./src
OBJDIR = $(SRCDIR)/obj
BINDIR = ./bin
LIBDIR = ./lib
INCDIR = ./include
TSTDIR = ./test

LIBS = -lgfortran -lm -lgomp
LIBFILE = -L$(LIBDIR) -lnblists

obj = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(1)))
src = $(addprefix $(SRCDIR)/, $(addsuffix .f90, $(1)))

OBJECTS = $(call obj,neighborLists lists global)

.PHONY: all clean install uninstall lib include test

.DEFAULT_GOAL := all

# Phony targets:

all: lib include

clean:
	rm -rf $(OBJDIR) $(LIBDIR) $(BINDIR) $(INCDIR)

install:
	cp $(LIBDIR)/libnblists.* $(PREFIX)/lib/
	cp $(INCDIR)/*.* $(PREFIX)/include/
	ldconfig

uninstall:
	rm -f $(PREFIX)/lib/libnblists.so
	rm -f $(addprefix $(PREFIX)/include/,nblists.h)
	ldconfig

lib: $(LIBDIR)/libnblists.so

include: $(INCDIR)/nblists.h $(INCDIR)/nblists.f90

test: $(addprefix $(TSTDIR)/,test_c_md)
	cd $(TSTDIR) && bash run_tests.sh

# Testing executables:

$(TSTDIR)/test_c_md: $(TSTDIR)/test_c_md.c $(INCDIR)/nblists.h $(LIBDIR)/libnblists.so
	$(CC) $(C_OPTS) -o $@ $(LN_OPTS) $< $(LIBFILE) -lm

$(TSTDIR)/%: $(TSTDIR)/%.c $(INCDIR)/nblists.h $(LIBDIR)/nblists.so
	$(FORT) $(F_OPTS) -o $@ $(LN_OPTS) -J$(OBJDIR) $< $(LIBFILE)

# Shared library and includes:

$(LIBDIR)/libnblists.so: $(OBJECTS)
	mkdir -p $(LIBDIR)
	$(FORT) $(F_OPTS) -shared -fPIC -o $@ $(OBJECTS) $(LIBS)

$(INCDIR)/nblists.h: $(SRCDIR)/nblists_header.h
	mkdir -p $(INCDIR)
	cp $< $@

$(INCDIR)/nblists.f90: $(SRCDIR)/nblists_module.f90
	mkdir -p $(INCDIR)
	cp $< $@

# Object files:

$(OBJDIR)/neighborLists.o: $(call src,neighborLists) $(call obj,lists global)
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/lists.o: $(SRCDIR)/lists.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/global.o: $(SRCDIR)/global.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<
