CC=			gcc
CFLAGS=		-Wall #-fno-inline-functions -fno-inline-functions-called-once
LIBS=		-lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all: CFLAGS+=-g -O3 -DNDEBUG
all: mfmi

debug: CFLAGS+=-DDEBUG -g -O0
debug: mfmi

mfmi:rle.o rope.o rlcsa.o main.o
		$(CC) $(CFLAGS) $(DFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -rf *.o mfmi
