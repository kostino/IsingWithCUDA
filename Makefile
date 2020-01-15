CC = gcc
NVCC = nvcc
CFLAGS = -O3 -Wall -Iinc
NVFLAGS = -O3 -Iinc 
TYPES = ising_V0 ising_V1 ising_V2 ising_V3

.PHONY: all lib clean

all: lib

lib: $(addsuffix .a, $(addprefix lib/ising_, $(TYPES)))

lib/%.a: obj/%.o
    mkdir -p lib
    ar rcs $@ $<

obj/ising_V0.o: src/ising_V0.c
    mkdir -p obj
    $(CC) $(CFLAGS) -o $@ -c $<

obj/ising_V1.o: src/ising_V1.cu
    mkdir -p obj
    $(NVCC) $(NVFLAGS) -o $@ -c $<

obj/ising_V2.o: src/ising_V2.cu
    mkdir -p obj
    $(NVCC) $(NVFLAGS) -o $@ -c $<

obj/ising_V3.o: src/ising_V3.cu
    mkdir -p obj
    $(NVCC) $(NVFLAGS) -o $@ -c $<

clean:
    rm obj/ising_*.o lib/ising_*.a
