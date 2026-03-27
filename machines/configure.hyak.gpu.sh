module use /sw/nvhpc/2024_245/modulefiles
module load cuda/12.4.1 ompi/4.1.4 cmake nvhpc-nompi/24.5

# Hyak-specific copy of the local Linux GPU configuration.
: "${PREFIX:=/mmfs1/gscratch/stf/abeck502/gkylsoft}"
./configure CC=nvcc CUDA_ARCH=70 --prefix=$PREFIX --cudamath-libdir=$LD_LIBRARY_PATH --use-mpi=yes --mpi-inc=/sw/ompi/4.1.4/include --mpi-lib=/sw/ompi/4.1.4/lib --use-nccl=yes --nccl-inc=/sw/nvhpc/2024_245/Linux_x86_64/24.5/comm_libs/nccl/include --nccl-lib=$LD_LIBRARY_PATH --use-cudss=yes --use-lua=yes
