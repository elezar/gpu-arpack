
FORT = gfortran -fPIC
CC = gcc -fPIC
OBJ_EXT = .o


LIB_ARPACK = -larpack
LIB_BLAS_LAPACK = -lblas -llapack
INC =

LIBS = $(LIB_ARPACK) $(LIB_BLAS_LAPACK)


MAGMA_PATH = /home/elezar/src/cuda/linalg/magma_1.1.0
CUDA_LIBS = -L/usr/local/cuda/lib64 -lcuda -lcublas -lcudart -L$(MAGMA_PATH)/lib -lmagmacontrol -lmagma -lmagmablas
MAGMA_INCLUDE = -I$(MAGMA_PATH)/include
NVCC_OPTIONS =

OUTPUT_DIR = ../lib


all: output_dir cpu_single

output_dir:
	mkdir -p $(OUTPUT_DIR)

arpack_drivers$(OBJ_EXT): arpack_drivers.f90
	$(FORT) -c arpack_drivers.f90

dense_cpu$(OBJ_EXT): dense_cpu.c common_functions.h interface.h
	$(CC) -c dense_cpu.c $(INC)

cpu_single: dense_cpu$(OBJ_EXT) arpack_drivers$(OBJ_EXT)
	$(CC) --shared -fPIC -fopenmp -o $(OUTPUT_DIR)/libsingle_dense_cpu.so $^ $(LIBS)

clean:
	rm -rf *$(OBJ_EXT)