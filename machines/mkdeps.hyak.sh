module load gcc/15.2.0 ompi/4.1.6
cd install-deps

# Hyak-specific dependency build helper.
: "${PREFIX:=/mmfs1/gscratch/stf/$USER/gkylsoft}"
./mkdeps.sh --build-openblas=yes --build-superlu=yes --build-luajit=yes --build-cudss=yes --prefix=$PREFIX
