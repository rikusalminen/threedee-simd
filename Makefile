CC=gcc
CFLAGS=-std=gnu99 -MMD -W -Wall -Wextra -pedantic
#CFLAGS+=-g -ggdb
CFLAGS+=-Iinclude/
CFLAGS+=-O3 -ffast-math
CFLAGS+=-march=native
ifeq ($(CC),gcc)
	CFLAGS+=-mfpmath=sse
endif
#CFLAGS+=-msse -msse2 -msse3 -msse4 -msse4.1 -msse4.2 -msse4a -msse2avx -mavx -mfma4
LDFLAGS=
LDLIBS=-lm
SRCS=$(wildcard *.c)
OBJS=$(SRCS:.c=.o)
DEPS=$(SRCS:.c=.d)

.PHONY: all clean
all: threedee-test

threedee-test: threedee-test.o

$(EXECUTABLE): $(OBJS)

-include $(DEPS)

clean:
	rm -f threedee-test
	rm -f $(OBJS) $(DEPS)

