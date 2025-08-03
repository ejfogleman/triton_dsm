# Makefile for the project

# Compiler and flags
CC = gcc
CFLAGS = -Wall -g
LDLIBS = -lm
TARGET = dsm
SRCS = dsm.c rngs.c

# Flags
CFLAGS_DEBUG   := -g -Wall
CFLAGS_RELEASE := -O3 -march=native -flto -DNDEBUG -Wall # -ffast-math
LDFLAGS := -lm

# Default target (debug)
all: debug

# Debug build
debug: CFLAGS := $(CFLAGS_DEBUG)
debug: $(TARGET)

# Release build
release: CFLAGS := $(CFLAGS_RELEASE)
release: $(TARGET)

# Build rule
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o $(TARGET) $(LDFLAGS)

# Clean up
.PHONY: all debug release clean
clean:
	rm -f $(TARGET)