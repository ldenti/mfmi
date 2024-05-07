CXX=		g++
CXXFLAGS=	-O3 -DNDEBUG -Wall -I./rlcsa/ -L./rlcsa/ -DMASSIVE_DATA_RLCSA -DMULTITHREAD_SUPPORT -D_GLIBCXX_PARALLEL -fopenmp
LIBS=		-lz -lrlcsa

.SUFFIXES:.c .o

.cpp.o:
		$(CXX) -c $(CXXFLAGS) $< -o $@

all: rl

rl:main_index.o main_exact.o main_pp.o main.o
		$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
		rm -rf *.o rl
