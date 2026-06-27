CC ?= gcc

ROOT := ../..
BUILD := build

CFLAGS := -std=c11 -Wall -Wextra -Wshadow -Wpedantic -Wconversion -Wsign-conversion -O1 -g
CFLAGS += -I$(ROOT) -I$(ROOT)/NIST

PROD_SRCS := \
  $(ROOT)/field.c \
  $(ROOT)/matrixmod.c \
  $(ROOT)/matrixelim.c \
  $(ROOT)/triform.c

TESTS := \
  test_layout \
  test_derived \
  test_corank \
  test_action \
  test_random

TARGETS := $(TESTS:%=$(BUILD)/%)

.PHONY: all run clean

all: $(TARGETS)

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/%: %.c test_common.h $(PROD_SRCS) | $(BUILD)
	$(CC) $(CFLAGS) $< $(PROD_SRCS) -o $@

run: all
	@for test in $(TARGETS); do \
	  echo "$$test"; \
	  ./$$test; \
	done

clean:
	rm -rf $(BUILD)
