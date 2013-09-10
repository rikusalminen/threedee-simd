# example with another compiler: make clean && make CC="/usr/local/bin/clang" && ./test

CC=gcc
CFLAGS=-std=gnu99 -MMD -W -Wall -Wextra -pedantic
CFLAGS+=-Wcast-align -Wcast-qual \
	-Wpointer-arith -Waggregate-return -Wunreachable-code \
	-Wfloat-equal -Wformat=2 -Wredundant-decls \
	-Wundef -Wdisabled-optimization -Wshadow \
	-Wmissing-braces -Wstrict-aliasing=2 -Wstrict-overflow=5 \
	-Wconversion -Wno-unused-parameter \
	-Wno-missing-field-initializers -Wno-missing-braces \
	-Wno-unused-function
CFLAGS+=-I../deps/
#CFLAGS+=-g -ggdb
CFLAGS+=-I../include/
CFLAGS+=-O3 -ffast-math
CFLAGS+=-march=native
ifeq ($(CC),gcc)
	CFLAGS+=-mfpmath=sse
endif
LDFLAGS=
LDLIBS=-lm
SRCS=$(wildcard *.c)
OBJS=$(SRCS:.c=.o)
DEPS=$(SRCS:.c=.d)

all: test

test: test.o

$(EXECUTABLE): $(OBJS)

-include $(DEPS)

clean:
	rm -f test
	rm -f $(OBJS) $(DEPS)

.PHONY: all clean
