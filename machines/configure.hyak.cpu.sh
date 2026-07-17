module purge
module load aocc/4.2.0 aocl/4.2.0 ompi/4.1.6-aocc openblas/0.3.30-aocc plasmawise/superlu/7.0.1-aocc plasmawise/luajit/2.1-aocc

# Hyak-specific copy of the local Linux CPU configuration. Installs to stf user directory by default (not home, but can make symbolic link there)
: "${PREFIX:=$HOME/gkylsoft}"

mkdir -p $PREFIX
ln -s /sw/ompi/4.1.6-4 $PREFIX/openmpi
ln -s /sw/openblas/0.3.30 $PREFIX/OpenBLAS
ln -s /sw/contrib/plasmawise-src/superlu-7.0.1 $PREFIX/superlu
ln -s /sw/contrib/plasmawise-src/luajit/2.1 $PREFIX/luajit

#./configure CC=clang --prefix=$PREFIX --use-mpi=yes --mpi-inc=/sw/ompi/4.1.6-4/include --mpi-lib=/sw/ompi/4.1.6-4/lib --use-lua=no --lua-inc=/sw/contrib/plasmawise-src/luajit/2.1/include --lua-lib=/sw/contrib/plasmawise-src/luajit/2.1/lib --lapack-inc=/sw/openblas/0.3.30/include --lapack-lib=/sw/openblas/0.3.30/lib --superlu-inc=/sw/contrib/plasmawise-src/superlu-7.0.1/include --superlu-lib=/sw/contrib/plasmawise-src/superlu-7.0.1/lib --app=moments

./configure CC=clang --prefix=$PREFIX --use-mpi=yes --use-lua=yes --app=moments
