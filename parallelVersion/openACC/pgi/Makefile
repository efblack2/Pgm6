# make file para la creacion de ObjectScan
#Name of the executable
PROJ = p6
EXT =
# end of Name of the executable

# Uncomment next line to build a double precision version
#DOUBLE=-DDOUBLE

############################################################
#  flags for the PGI pgcc compiler
    CC=pgcc
	#OMP= -mp
	ACC= -acc -ta=tesla:fastmath,cc50,lineinfo -Mcuda=ptxinfo #-ta=nvidia,keepgpu,keepptx,nollvm
	COMPOP= -fast -c11  -Minfo=all  -Mvect=simd $(OMP) $(ACC) $(DOUBLE)
	#ACC= -acc -ta=tesla:fastmath,cc30,cuda7.5,lineinfo  -Mcuda=ptxinfo,maxregcount:25 #-ta=nvidia,keepgpu,keepptx,nollvm
	#COMPOP= -O3 -fast -c11  -Minfo=all  -Mvect=simd $(OMP) $(ACC)
    # -mp -Mvect=simd or -Mnovect
	PROFILE= -Mprof=ccff #-Mprof=lines
	COMPILER_FLAGS= -c $(COMPOP)  $(PROFILE)
	BIND_FLAGS=  $(PROFILE)
############################################################

info :
	@echo "The folowing are possible:"
	@echo " "make $(PROJ)$(EXT) .... to build the program named $(PROJ)$(EXT)
	@echo " make clean"


.SUFFIXES: .o .c

.PHONY : clean

all: $(PROJ)$(EXT)

############# Changes are done here #############################################
objs= \
    ./main.o \
    ./advection.o \
    ./advec1.o \
    ./pgf.o \
    ./diffusion.o \
    ./bc.o \
    ./bc4T.o \
    ./ic.o \
    ./update.o \
    ./getDensity.o \
	./putfield.o

#    ./getTakacsError.o \

# place to put the dependencies


./main.o: \
    ./main.c\
    ./dimCube.h \
    ./real.h \
    ./prototypes.h

./advection.o: \
    ./advection.c \
    ./real.h

./advec1.o: \
    ./advec1.c \
    ./real.h

./diffusion.o: \
    ./diffusion.c \
    ./real.h

./pgf.o: \
    ./pgf.c \
    ./real.h


./bc.o: \
    ./bc.c \
    ./real.h

./bc4T.o: \
    ./bc4T.c \
    ./real.h

./ic.o: \
    ./ic.c \
    ./real.h



./putfield.o: \
    ./putfield.c \
    ./real.h


./update.o: \
    ./update.c \
    ./real.h

./getDensity.o: \
    ./getDensity.c \
    ./real.h


#./getTakacsError.o: \
#    ./getTakacsError.c




#######################  do not change anything from here on ####################


%.o : %.c
	$(CC) $(COMPILER_FLAGS) $<

$(PROJ)$(EXT): $(objs)
	$(CC) $(objs)  -o  $@ $(BIND_FLAGS) $(COMPOP)

clean:
	@rm -rf $(objs)  *.o

