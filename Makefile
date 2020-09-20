GCC=/usr/bin/gcc
CFLAGS=--std=c99 -Wall -Wextra -Werror -pedantic-errors
DEBUG_FLAGS=-O0 -g -D__DEBUG__ -pg
SOURCES=$(wildcard *.c)
OBJECTS=$(SOURCES:.c=.o)
LIBS=m
LIBFLAGS=$(addprefix -l, $(LIBS))
EXEC=cluster

.PHONY: all clean test debug
all: $(EXEC)

debug: CFLAGS += $(DEBUG_FLAGS)
debug: $(EXEC)

$(EXEC): $(OBJECTS)
	$(GCC) $(CFLAGS) $^ -o $@ $(LIBFLAGS)

%.o: %.c
	$(GCC) -c $(CFLAGS) $^ -o $@

clean:
	rm -f $(OBJECTS) $(EXEC)
