CC=			gcc
CFLAGS=		-Wall -O0 -g #-fno-inline-functions -fno-inline-functions-called-once
DFLAGS=
INCLUDES=	
LIBS=		-lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:ppp

ppp:rle.o rlcsa.o main.o
		$(CC) $(CFLAGS) $(DFLAGS) $^ -o $@ $(LIBS)

rle.o:rle.h
rlcsa.o:rle.h rlcsa.h 
main.o:rle.h rlcsa.h

clean:
		rm -rf *.o ppp
