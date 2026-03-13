# Bio-CFD

"Description of what the code does"

## Prerequisites

To build this project, we recommend using [`fpm`][1].  Instructions
for installing `fpm` can be found [here][2].

## Building

Once `fpm` has been installed, building the code becomes
simple. First, decide on your compiler and compiler options by setting
the `FPM_FC` and `FPM_FFLAGS` environment variables e.g. to use
`nvfortran` and target `cuda` offloading do

```
export FPM_FC="nvfortran"
export FPM_FFLAGS="-Mr8 -acc -mcmodel=medium -gpu=managed -cuda -mp"
```

You can then simply run

```
fpm build
```
in the root directory of the project.

## Running the code

Running the code is also simple. Simply execute the command

```
fpm run
```

inside the root directory of the project. Be sure to set the
environment variables described in the build section first. `fpm` is
smart enough to know if you have changed the code since building, and
will automatically call `fpm build` before `fpm run` for you if
needed.


[1]: https://fpm.fortran-lang.org/index.html
[2]: https://fpm.fortran-lang.org/install/index.html#install

## Acknowledgements

This software was substantially developed with research software engineering support funded by the Baskerville Tier 2 HPC service (https://www.baskerville.ac.uk/). Baskerville was funded by the EPSRC and UKRI through the World Class Labs scheme (EP/T022221/1) and the Digital Research Infrastructure programme (EP/W032244/1) and is operated by Advanced Research Computing at the University of Birmingham.
