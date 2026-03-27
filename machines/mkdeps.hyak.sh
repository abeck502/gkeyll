module load ompi/4.1.4 cmake
cd install-deps

# Hyak-specific dependency build helper.
: "${PREFIX:=/mmfs1/gscratch/stf/abeck502/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --build-luajit=yes --build-cudss=yes --prefix=$PREFIX
