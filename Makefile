GCC=/usr/bin/gcc
CFLAGS=-ansi -Wall -Wextra -Werror -pedantic-errors
DEBUG_FLAGS=-O0 -g
SOURCES=$(wildcard *.c)
# SOURCES=$(filter-out adjacency_matrix.c, $(wildcard *.c))
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
