
FORT = gfortran -fPIC
CC = gcc -fPIC
OBJ_EXT=.o


LIB_ARPACK=-larpack
ACML_PATH=/opt/acml5.1.0/gfortran64
INC=-I$(ACML_PATH)/include
LIB_BLAS_LAPACK=-L$(ACML_PATH)/lib -lacml

LIBS=$(LIB_ARPACK) $(LIB_BLAS_LAPACK)
#LIB_BLAS_LAPACK=-lblas -llapack
MAGMA_PATH=/home/elezar/src/cuda/linalg/magma_1.1.0
CUDA_LIBS=-L/usr/local/cuda/lib64 -lcuda -lcublas -lcudart -L$(MAGMA_PATH)/lib -lmagmacontrol -lmagma -lmagmablas
MAGMA_INCLUDE =-I$(MAGMA_PATH)/include
NVCC_OPTIONS=-ccbin /home/elezar/src/cuda/utils/gcc44

LIB_DIR=../lib




all: cpu_single

arpack_drivers$(OBJ_EXT): arpack_drivers.f90
	$(FORT) -c arpack_drivers.f90

dense_cpu$(OBJ_EXT): dense_cpu.c common_functions.h interface.h
	$(CC) -c dense_cpu.c $(INC)

cpu_single: dense_cpu$(OBJ_EXT) arpack_drivers$(OBJ_EXT)
	$(CC) --shared -fPIC -fopenmp -o $(LIB_DIR)/libsingle_dense_cpu.so $^ $(LIBS)

clean:
	rm -rf *$(OBJ_EXT)