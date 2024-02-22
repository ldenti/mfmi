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

rle.o:rle.h
rope.o:rle.h rope.h
rlcsa.o:rope.h kvec.h rlcsa.h
main.o:rlcsa.h

clean:
		rm -rf *.o mfmi
