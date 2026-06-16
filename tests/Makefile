SHELL := /bin/bash

CC := gcc
LIBS = -lssl -lcrypto

ifdef DEBUG
CFLAGS := -g -Wall -DDEBUG
OBJDIR := build/debug
else
CFLAGS := -O3 -Wall
OBJDIR := build
endif

CFLAGS += -I.. -I../NIST

CORE_OBJECTS = \
	$(OBJDIR)/meds.o \
	$(OBJDIR)/util.o \
	$(OBJDIR)/seed.o \
	$(OBJDIR)/osfreq.o \
	$(OBJDIR)/fips202.o \
	$(OBJDIR)/matrixmod.o \
	$(OBJDIR)/bitstream.o \
	$(OBJDIR)/randombytes.o

HEADERS = $(wildcard ../*.h)

.PHONY: all test clean

all: test

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: ../%.c $(HEADERS) | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/test: test.c $(CORE_OBJECTS)
	$(CC) test.c $(CORE_OBJECTS) $(CFLAGS) $(LIBS) -o $@

test: $(OBJDIR)/test

clean:
	rm -rf build
