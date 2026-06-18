module purge
module load gcc/15.2.0 ompi/4.1.6 

# Hyak-specific copy of the local Linux CPU configuration.
# Keep machines/configure.linux.cpu.sh aligned with upstream defaults.
: "${PREFIX:=/mmfs1/gscratch/stf/abeck502/gkylsoft}"
./configure CC=cc --prefix=$PREFIX --use-mpi=yes --mpi-inc=/sw/ompi/4.1.6/include --mpi-lib=/sw/ompi/4.1.6/lib --use-lua=yes --app=moments
