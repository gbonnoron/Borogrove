
CC=g++
CFLAGS=-c -Wall -std=c++11 -Ofast -march=native -mtune=native -fno-schedule-insns -funroll-loops -ffinite-math-only -g -ggdb
LDFLAGS=-march=native -mtune=native -lfftw3
SOURCES= ../circulant_ring.cpp \
		../rlwe.cpp \
		../operations.cpp \
		../ksw_key.cpp \
		../rgsw.cpp \
		../fft.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: test_lwe test_gsw test_ops test_int

check: all;
	./test_lwe && ./test_gsw && ./test_ops

test_lwe: $(OBJECTS) test_lwe.o
	$(CC) $(OBJECTS) test_lwe.o -o $@ $(LDFLAGS)

test_gsw: $(OBJECTS) test_gsw.o
	$(CC) $(OBJECTS) test_gsw.o -o $@ $(LDFLAGS)

test_ops: $(OBJECTS) test_ops.o
	$(CC) $(OBJECTS) test_ops.o -o $@ $(LDFLAGS)

test_int: $(OBJECTS) test_int.o
	$(CC) $(OBJECTS) test_int.o -o $@ $(LDFLAGS)

test_mult:
	$(CC) test_mult.cpp -o $@ $(CFLAGS) $(LDFLAGS)

clean:
	rm -f test_lwe test_gsw test_ops test_int *.o

rebuild: clean;
	make -j5

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@
