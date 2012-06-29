CFLAGS=-std=gnu99 -Iinclude/ -W -Wall -Wextra -pedantic -O3 -march=native
#CFLAGS+=-msse -msse2 -msse3 -msse4 -msse4.1 -msse4.2 -msse4a -msse2avx -mavx -mfma4
LIBS=-lm
CC=gcc
SRCS=$(wildcard *.c)
OBJS=$(SRCS:.c=.o)
DEPS=$(SRCS:.c=.d)
EXECUTABLE=threedee-test

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) $(LIBS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.d: %.c
	$(CC) $(CFLAGS) -MM -o $@ $<

-include $(DEPS)

.PHONY: clean

clean:
	rm -f $(EXECUTABLE)
	rm -f $(OBJS) $(DEPS)

