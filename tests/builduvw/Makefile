ifeq ($(origin CC), default)
CC := gcc
endif

ROOT := ../..
BUILD := build

BASE_CFLAGS := -std=c11 -O1 -g
WARN_CFLAGS := \
  -Wall \
  -Wextra \
  -Wshadow \
  -Wpedantic \
  -Wconversion \
  -Wsign-conversion
INCLUDES := -I$(ROOT) -I$(ROOT)/NIST
SAN_CFLAGS := -fsanitize=address,undefined -fno-omit-frame-pointer

CFLAGS ?= $(BASE_CFLAGS) $(WARN_CFLAGS) $(INCLUDES)
PROD_CFLAGS ?= $(BASE_CFLAGS) -Wall -Wextra -Wshadow -Wpedantic $(INCLUDES)
LDFLAGS ?=

PROD_OBJECTS := \
  $(BUILD)/field.o \
  $(BUILD)/matrixmod.o \
  $(BUILD)/matrixelim.o \
  $(BUILD)/triform.o \
  $(BUILD)/canonical.o

TESTS := \
  test_input \
  test_success \
  test_failure \
  test_projective \
  test_random

TARGETS := $(TESTS:%=$(BUILD)/%)

.PHONY: all run strict sanitize clean

all: $(TARGETS)

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/field.o: $(ROOT)/field.c | $(BUILD)
	$(CC) $(PROD_CFLAGS) -c $< -o $@

$(BUILD)/matrixmod.o: $(ROOT)/matrixmod.c | $(BUILD)
	$(CC) $(PROD_CFLAGS) -c $< -o $@

$(BUILD)/matrixelim.o: $(ROOT)/matrixelim.c | $(BUILD)
	$(CC) $(PROD_CFLAGS) -c $< -o $@

$(BUILD)/triform.o: $(ROOT)/triform.c | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD)/canonical.o: $(ROOT)/canonical.c | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD)/%.test.o: %.c test_common.h test_reference.h | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD)/%: $(BUILD)/%.test.o $(PROD_OBJECTS)
	$(CC) $^ $(CFLAGS) $(LDFLAGS) -o $@

run: all
	@for test in $(TARGETS); do \
	  echo "$$test"; \
	  ./$$test; \
	done

strict:
	$(MAKE) clean
	$(MAKE) \
	  CFLAGS="$(BASE_CFLAGS) $(WARN_CFLAGS) $(INCLUDES) -Werror" \
	  PROD_CFLAGS="$(BASE_CFLAGS) -Wall -Wextra -Wshadow -Wpedantic $(INCLUDES)" \
	  run

sanitize:
	$(MAKE) clean
	ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 \
	UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
	$(MAKE) \
	  CFLAGS="$(BASE_CFLAGS) $(WARN_CFLAGS) $(INCLUDES) $(SAN_CFLAGS)" \
	  PROD_CFLAGS="$(BASE_CFLAGS) -Wall -Wextra -Wshadow -Wpedantic $(INCLUDES) $(SAN_CFLAGS)" \
	  run

clean:
	rm -rf $(BUILD)
