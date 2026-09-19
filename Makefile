# Makefile
CC      = gcc
# -include at the end of this file can supply a default goal; pin it.
.DEFAULT_GOAL := all
# -MMD -MP: per-object header deps, so a header change rebuilds its users.
CFLAGS = -Wall -Wextra -O2 -Isrc -Ithird_party -MMD -MP
ifdef DEBUG
CFLAGS = -Wall -Wextra -g -O0 -Isrc -Ithird_party -MMD -MP
endif
LIBS    = -lm

# Dataset path
DATA_IMU = ~/datasets/euroc/MH_01_easy/mav0/imu0/data.csv
DATA_GT  = ~/datasets/euroc/MH_01_easy/mav0/state_groundtruth_estimate0/data.csv
DATA_CAM = ~/datasets/euroc/MH_01_easy/mav0/cam0

SRCS    = $(wildcard src/*.c)
OBJS    = $(patsubst src/%.c,build/%.o,$(SRCS))

# Module objects only (no main.o) — linked into every test
MOD_OBJS  = $(filter-out build/main.o,$(OBJS))

# Tests: every tests/*.c becomes its own binary build/test_*
TEST_SRCS = $(wildcard tests/*.c)
TEST_BINS = $(patsubst tests/%.c,build/%,$(TEST_SRCS))

# Header deps from -MMD. Must be included at the END of this file.
DEPS = $(OBJS:.o=.d) $(TEST_BINS:=.d)

all: build/main

build:
	mkdir -p build

build/main: $(OBJS)
	$(CC) $(OBJS) -o $@ $(LIBS)

build/%.o: src/%.c | build
	$(CC) $(CFLAGS) -c $< -o $@

build/%: tests/%.c $(MOD_OBJS) | build
	$(CC) $(CFLAGS) $< $(MOD_OBJS) -o $@ $(LIBS)

test: $(TEST_BINS)
	@for t in $(TEST_BINS); do echo "== $$t"; $$t || exit 1; done

run: build/main
	@./build/main $(DATA_IMU) $(DATA_GT) $(DATA_CAM)

clean:
	rm -rf build

.PHONY: all test run clean

-include $(DEPS)

