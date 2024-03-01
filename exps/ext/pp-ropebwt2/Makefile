CC=		gcc
CFLAGS=		-O3 -DNDEBUG -Wall -I./ropebwt2/ -fopenmp
LIBS=		-lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $< -o $@

all: rt2

rt2:ropebwt2/rld0.o ropebwt2/mrope.o ropebwt2/rope.o ropebwt2/rle.o main_index.o main_pp.o main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -rf *.o rt2
