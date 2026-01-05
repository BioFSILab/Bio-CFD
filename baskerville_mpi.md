## MPI on Baskerville

To build and run using MPI the process is currently a bit more
complicated than the non-MPI version.

### Installing

Rather than using `fpm build` we will use `fpm install` to install the
MPI-enabled executable into a path of your choice. 

As a starting point, we will use the same settings we used for the
non-MPI enabled version of the code i.e. run
```
export FPM_FC="nvfortran"
export FPM_FFLAGS="-Mr8 -acc -mcmodel=medium -gpu=managed -cuda -mp"
```
we now need to add the appropiate MPI flags for compilation. Run
```
MPI_FFLAGS=$(mpif90 -showme | cut -d ' ' -f 2-)
```
and then finally 
```
export FPM_FFLAGS="${FPM_FFLAGS} ${MPI_FFLAGS} -DBIOCFD_MPI"
```
If you now look at `FPM_FFLAGS` with e.g.
```
echo ${FPM_FFLAGS}
```
you should see the original flags we set, followed by several paths
and library with MPI in the name, and finally `-DBIOCFD_MPI`.

At this stage, we want to load all the required modules using the
module system.

On Baskerville these are
```
module load baskerville
module load HDF5/1.14.0-NVHPC-25.7-CUDA-12.6.0-serial
module load nvompi/2025.07
```

Once this is done we can run `fpm install` with a path to where we
want to build the executable. For example,
```
fpm install --profile release --prefix mpi_build/
```

A complete SLURM script to build the executable might look like (be
sure to change the two paths marked with `<>` to point to the correct
places)

```
#!/bin/bash
#SBATCH --qos=bham
#SBATCH --account=bosec-butterflies
#SBATCH --nodes 1
#SBATCH --time 00:10:00
#SBATCH --cpus 18

module purge
module load baskerville
module load HDF5/1.14.0-NVHPC-25.7-CUDA-12.6.0-serial
module load nvompi

# TODO: Replace with the path to the FPM installation if not already in your $PATH
export PATH="<PATH_TO_WHERE_FPM_IS_INSTALLED>:${PATH}"

export FPM_FC="nvfortran"
export FPM_FFLAGS="-Mr8 -acc -mcmodel=medium -gpu=mem:managed -cuda -Minfo=accel,inline -mp -O3 -gpu=cc80,cc90 -tp=x86-64-v3"
MPI_FFLAGS=$(mpif90 -showme | cut -d ' ' -f 2-)
export FPM_FFLAGS="${FPM_FFLAGS} ${MPI_FFLAGS} -DBIOCFD_MPI"

# TODO: Replace with the path where you want to install the executable
fpm install --profile release --prefix <PATH_TO_PLACE_THE_EXECUTABLE>
```

After running this you will find an executable in `mpi_build/bin`
called `main`. Make a note of this path, because we will put this into
our SLURM submission file.

### Running

Navigate to the directory where you want to run the program. This will
be where you have put all the input files and an empty `out`
directory.

In this directory create a new script as below

```
#!/bin/bash

#SBATCH --qos=bham
#SBATCH --account=bosec-butterflies
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 36

module purge
module load baskerville
module load HDF5/1.14.0-NVHPC-25.7-CUDA-12.6.0-serial
module load nvompi

PROG=<FULL_PATH_TO_THE_EXECUTABLE_YOU_CREATED_IN_THE_PREVIOUS_STEP>
srun -o "slurm-%j-%N.out" ${PROG}
```
If we call this `submit.sh` then we can submit this with
```
sbatch --gpus N --time HH:MM:SS --constraint a100_80 submit.sh
```
Replace `N` with the number of GPUs you want to use and `HH:MM:SS`
with the time in hours, minutes, and seconds, that you'll need for
your job. `--constraint a100_80` limits us to only 80GB A100 GPUs, as
large grids can use up a lot of memory.
