# Borogrove
Large FHE gates using Homomorphic Accumulator

## How to cite ##

	@unpublished{he8,
	    Note = {Available at \url{https://github.com/gbonnoron/Borogrove}},
	    Title = {Large FHE Gates using Homomorphic Accumulator},
	    Author = {Bonnoron, Guillaume \and Ducas, Léo \and Fillinger, Max},
	    Year = 2017,
	    url = {https://github.com/gbonnoron/Borogrove}
	}

# Compilation #

## Dependencies ##

### Required ###

- FFTW 3.3.6-pl2, or higher [http://www.fftw.org/](http://www.fftw.org/)
- g++ 

## Installation ##

You should download the source code from github and then run

    make

By default, this tells FFTW to use sub-optimal target (`FFTW_MEASURE`) when doing the precomputations.
You may change this to optimal target (`FFTW_PATIENT`) by editing the first lines of the Makefile.

    #PLAN_MODE=FFTW_ESTIMATE # worst performance, fastest to pre-compute
    PLAN_MODE=FFTW_MEASURE
    #PLAN_MODE=FFTW_PATIENT # best performance, slowest to pre-compute

Leave the mode you want, and only this one, uncommented. Please refer to the FFTW manual for further details about plans and wisdoms.
Use `make rebuild` if you changed the plan mode and want to recompile the whole project

## Check ##

Type

	make check


## Optimization ##

The default compilation flag are `-Ofast -march=native -mtune=native -fno-schedule-insns -funroll-loops -ffinite-math-only`. You may change them.


# How to use #

## Main program ##

To start the main program, run

    ./he8

You may provide a command line argument to specify the number of gates you want to evaluate.
Note: At the first launch, FFTW will generate wisdom, which can take some time depending of the selected mode.

## Tests and stats ##
In case you want to modify parameters and check correctness of the scheme, you can run the different unit tests in the `./tests` folder.

    cd tests/
    make

Be sure you duplicate the `wisdom*` files you may have in the main directory to avoid doing the whole fftw precomputation again.

In the `stats` folder, you have one program for each elementary operations and one global program that covers the whole scheme. These programs execute the operation(s) and measure the output error.
The input error variances for each operation is set, by default, to the expected output value of the operation that comes just before during a complete gate execution, for the given parameter sets.

    cd stats/
    make
    ./global

Again, you can duplicate your wisdoms here to save time.
