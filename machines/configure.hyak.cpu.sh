module load cmake ompi/4.1.4

# Hyak-specific copy of the local Linux CPU configuration.
# Keep machines/configure.linux.cpu.sh aligned with upstream defaults.
: "${PREFIX:=/mmfs1/gscratch/stf/abeck502/gkylsoft}"
./configure CC=cc --prefix=$PREFIX --use-mpi=yes --mpi-inc=/sw/ompi/4.1.4/include --mpi-lib=/sw/ompi/4.1.4/lib --use-lua=yes --app=moments
