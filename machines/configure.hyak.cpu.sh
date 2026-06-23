module purge
module load gcc/15.2.0 ompi/4.1.6 

# Hyak-specific copy of the local Linux CPU configuration. Installs to stf user directory by default (not home, but can make symbolic link there)
: "${PREFIX:=/mmfs1/gscratch/stf/$USER/gkylsoft}"
./configure CC=cc --prefix=$PREFIX --use-mpi=yes --mpi-inc=/sw/ompi/4.1.6/include --mpi-lib=/sw/ompi/4.1.6/lib --use-lua=yes --app=moments
