SHELL := /bin/bash

LIBS = -lssl -lcrypto
CC := gcc

ifdef DEBUG
CFLAGS := -g -Wall -DDEBUG
OBJDIR := debug
else
CFLAGS := -O3 -Wall
OBJDIR := build
endif

ifdef LOG
CFLAGS += -DDEBUG
endif

CFLAGS += -I$(OBJDIR) -INIST

EXES = bench PQCgenKAT_sign
SPEED_EXE = $(OBJDIR)/test_speed

ifneq ($(wildcard test_meds2end_seedless.c),)
EXES += test_meds2end_seedless
endif

TARGETS := ${EXES:%=$(OBJDIR)/%}

.PHONY: default clean KAT BENCH SPEED

default: $(EXES)

OBJECTS = meds.o util.o osfreq.o fips202.o field.o matrixmod.o matrixelim.o triform.o canonical.o corank1.o trine_expand.o trine_codec.o bitstream.o randombytes.o
HEADERS = $(wildcard *.h)

BUILDOBJ := ${OBJECTS:%=$(OBJDIR)/%}


$(EXES) : % :
	@make $(OBJDIR)/$(@F)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BUILDOBJ) : $(OBJDIR)/%.o: %.c $(HEADERS) | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGETS) : $(OBJDIR)/%: %.c $(BUILDOBJ)
	$(CC) $(@F).c $(BUILDOBJ) $(CFLAGS) $(LIBS) -o $@

$(SPEED_EXE): test_speed/test_speed.c test_speed/cpucycles.c test_speed/speed_print.c $(BUILDOBJ)
	$(CC) test_speed/test_speed.c test_speed/cpucycles.c test_speed/speed_print.c $(BUILDOBJ) $(CFLAGS) -Itest_speed $(LIBS) -o $@

KAT: PQCgenKAT_sign
	$(OBJDIR)/PQCgenKAT_sign

BENCH: bench
	$(OBJDIR)/bench

SPEED: $(SPEED_EXE)
	$(SPEED_EXE)

clean:
	rm -rf build/ debug/ PQCsignKAT_*
