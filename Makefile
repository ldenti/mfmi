CC=			gcc
CFLAGS=		-Wall #-fno-inline-functions -fno-inline-functions-called-once
DFLAGS=
INCLUDES=	
LIBS=		-lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all: CFLAGS+=-g -O3 -DNDEBUG
all: ppp

debug: CFLAGS+=-DDEBUG -g -O0
debug: ppp

ppp:rle.o rlcsa.o main.o
		$(CC) $(CFLAGS) $(DFLAGS) $^ -o $@ $(LIBS)

# brle.o:brle.h
rle.o:rle.h
rlcsa.o:rle.h kvec.h rlcsa.h
main.o:rle.h rlcsa.h

clean:
		rm -rf *.o ppp
