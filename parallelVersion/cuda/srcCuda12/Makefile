# make file para la creacion de ObjectScan
#Name of the executable
PROJ = p6
EXT = 
# end of Name of the executable

# Uncomment next line to build a double precision version
#DOUBLE=-DDOUBLE

############################################################
#  #flags for the GNU g++ compiler
    CC=gcc
    # -fopenmp  -fno-tree-vectorize  -ftree-vectorize -ffast-math
    #PROFILE= -p -g  
#  flags for the nvcc compiler
    CC=gcc
    CUDA=nvcc
    LINKER=nvcc
    #PROFILE= -pg -lineinfo
    #  -maxrregcount=32 
    EXTRA_CUDACFLAGS= -std=c++11 -rdc=true -maxrregcount=32  --use_fast_math -lcudart --ptxas-options=-v -gencode arch=compute_50,code=sm_50
    COMPILER_FLAGS= -c -O3 -ccbin=gcc   $(DOUBLE) $(PROFILE)
    CUDACOMPILER_FLAGS= $(COMPILER_FLAGS)  $(EXTRA_CUDACFLAGS)
    BIND_FLAGS= -lm -lcudart $(PROFILE)    
############################################################

info :
	@echo "The folowing are possible:"
	@echo " "make $(PROJ)$(EXT) .... to build the program named $(PROJ)$(EXT)
	@echo " make clean"


.SUFFIXES: .cu .c .o

.PHONY : clean

all: $(PROJ)$(EXT)

############# Changes are done here #############################################
objs= \
    ./main.o \
    ./advectionTB.o \
    ./advectionNS.o \
    ./advectionEW.o \
    ./advecUVW.o \
    ./pgf1.o \
    ./pgf2.o \
    ./diffusion.o \
    ./bcNS.o \
    ./bcEW.o \
    ./bcTB.o \
    ./ic.o \
    ./update.o \
    ./getDensity.o \
	./putfield.o
     
 #    ./getTakacsError.o \

# place to put the dependencies


./main.o: \
    ./main.cu\
    ./dimCube.h \
    ./real.h \
    ./std2file.hpp \
    ./fatal.hpp \
    ./constant.cuh \
    ./prototypes.cuh

./advectionTB.o: \
    ./advectionTB.cu \
    ./constant.cuh \
    ./macros.h \
    ./real.h

./advectionNS.o: \
    ./advectionNS.cu \
    ./constant.cuh \
    ./macros.h \
    ./real.h

./advectionEW.o: \
    ./advectionEW.cu \
    ./constant.cuh \
    ./macros.h \
    ./real.h

./advecUVW.o: \
    ./advecUVW.cu \
    ./constant.cuh \
    ./macros.h \
    ./real.h


./diffusion.o: \
    ./diffusion.cu \
    ./constant.cuh \
    ./macros.h \
    ./real.h
    
./pgf1.o: \
    ./pgf1.cu \
    ./constant.cuh \
    ./macros.h \
    ./real.h

./pgf2.o: \
    ./pgf2.cu \
    ./constant.cuh \
    ./macros.h \
    ./real.h

./bcNS.o: \
    ./bcNS.cu \
    ./constant.cuh \
    ./macros.h \
    ./real.h

./bcEW.o: \
    ./bcEW.cu \
    ./constant.cuh \
    ./macros.h \
    ./real.h

./bcTB.o: \
    ./bcTB.cu \
    ./constant.cuh \
    ./macros.h \
    ./real.h

./ic.o: \
    ./ic.cu \
    ./real.h



./putfield.o: \
    ./putfield.cu \
    ./real.h


#./getTakacsError.o: \
#    ./getTakacsError.c

./update.o: \
    ./update.cu \
    ./real.h
    
./getDensity.o: \
    ./getDensity.cu \
    ./real.h
    

    

#######################  do not change anything from here on ####################

%.o : %.c
	$(CC)  $(COMPILER_FLAGS) -o $*.o  $<

%.o : %.cu
	$(CUDA)  $(CUDACOMPILER_FLAGS) -o $*.o  $<

$(PROJ)$(EXT): $(objs)
	$(LINKER) $(BIND_FLAGS) -o  $@ $(objs)

clean:
	@rm -rf $(objs)  *.o

