## Join output and error in a single file
#PBS -j oe
## Export the environment vaiables from your shell
#PBS -V
#PBS -l nodes=node
cd $PBS_O_WORKDIR
./aout
