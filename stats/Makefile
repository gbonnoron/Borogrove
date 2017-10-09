
CC=g++
CFLAGS=-g -ggdb -Wall -std=c++11 -O3 -march=native -mtune=native
LDFLAGS=-O3 -march=native -mtune=native -lfftw3
SOURCES= ../circulant_ring.cpp \
		../rlwe.cpp \
		../operations.cpp \
		../ksw_key.cpp \
		../rgsw.cpp \
		../fft.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all: 01-add 02-mod_switch 03-ext_exp_inner 04-mod_switch 05-exp_crt 06-mod_switch 07-fun_exp_extract 11-ext_mult 12-key_switch global

01-add: $(OBJECTS) 01-add.o
	$(CC) $(OBJECTS) 01-add.o -o $@ $(LDFLAGS)

02-mod_switch: $(OBJECTS) 02-mod_switch.o
	$(CC) $(OBJECTS) 02-mod_switch.o -o $@ $(LDFLAGS)

03-ext_exp_inner: $(OBJECTS) 03-ext_exp_inner.o
	$(CC) $(OBJECTS) 03-ext_exp_inner.o -o $@ $(LDFLAGS)

04-mod_switch: $(OBJECTS) 04-mod_switch.o
	$(CC) $(OBJECTS) 04-mod_switch.o -o $@ $(LDFLAGS)

05-exp_crt: $(OBJECTS) 05-exp_crt.o
	$(CC) $(OBJECTS) 05-exp_crt.o -o $@ $(LDFLAGS)

06-mod_switch: $(OBJECTS) 06-mod_switch.o
	$(CC) $(OBJECTS) 06-mod_switch.o -o $@ $(LDFLAGS)

07-fun_exp_extract: $(OBJECTS) 07-fun_exp_extract.o
	$(CC) $(OBJECTS) 07-fun_exp_extract.o -o $@ $(LDFLAGS)

11-ext_mult: $(OBJECTS) 11-ext_mult.o
	$(CC) $(OBJECTS) 11-ext_mult.o -o $@ $(LDFLAGS)

12-key_switch: $(OBJECTS) 12-key_switch.o
	$(CC) $(OBJECTS) 12-key_switch.o -o $@ $(LDFLAGS)

global: $(OBJECTS) global.o
	$(CC) $(OBJECTS) global.o -o $@ $(LDFLAGS)

rebuild: clean;
	make -j5

clean:
	rm -f *.o 01-add 02-mod_switch 03-ext_exp_inner 04-mod_switch 05-exp_crt 06-mod_switch 07-fun_exp_extract 11-ext_mult 12-key_switch global

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@
