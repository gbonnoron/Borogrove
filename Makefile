
#PLAN_MODE=FFTW_ESTIMATE # worst performance, fastest to pre-compute
#PLAN_MODE=FFTW_MEASURE
PLAN_MODE=FFTW_PATIENT # best performance, slowest to pre-compute

CC=g++
CFLAGS=-DPLAN_MODE=$(PLAN_MODE)
CFLAGS+=-c -Wall -std=c++11 -Ofast -march=native -mtune=native -fno-schedule-insns -funroll-loops -ffinite-math-only -g -ggdb
#CFLAGS=-c -Wall -std=c++11 -O3 -march=native -mtune=native -pg -g -ggdb
LDFLAGS=-march=native -mtune=native -lfftw3
SOURCES= he8.cpp \
		circulant_ring.cpp \
		rlwe.cpp \
		operations.cpp \
		ksw_key.cpp \
		rgsw.cpp \
		fft.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=he8

all: $(SOURCES) $(EXECUTABLE) 

check:
	make check -C tests/

tests:
	make -C tests/

stats:
	make -C stats/

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)

mr_proper:
	rm -f $(EXECUTABLE) $(OBJECTS)
	make clean -C tests/
	make clean -C stats/
	rm -fr doc

doc: $(SOURCES)
	/usr/bin/doxygen Doxyfile
	make -C doc/latex

rebuild: clean;
	make -j5

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
