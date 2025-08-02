# Makefile for the project

CC = gcc
CFLAGS = -Wall -g
LDLIBS = -lm
TARGET = dsm
SRCS = dsm.c rngs.c

# Phony targets

.PHONY: all, clean

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LDLIBS)

clean:
	rm -f $(TARGET)
