CC=			gcc
CXX=		g++
CFLAGS=		-Wall -fopenmp #-fno-inline-functions -fno-inline-functions-called-once
INCLUDES=	
CXXFLAGS=	$(CFLAGS)
LIBS=		-lz # -lpthread

VPATH = bits:misc

.c.o:
		$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

all: CFLAGS+=-g -O2 -DNDEBUG
all: CXXFLAGS+=-g -O2 -DNDEBUG
all: mfmi

debug: CFLAGS+=-DDEBUG -g -O0
debug: CXXFLAGS+=-DDEBUG -g -O0
debug: mfmi

mfmi:bitvector.o bitbuffer.o rlevector.o rle.o rope.o rlcsa.o rld0.o main_index.o main.o # main_search.o main_pp.o
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -rf *.o mfmi
