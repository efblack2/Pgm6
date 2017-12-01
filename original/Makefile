# make file para la creacion de ObjectScan
#Name of the executable
PROJ = p6
EXT = 
# end of Name of the executable

############################################################
#  flags for the GNU g++ compiler
   CC=gcc
	COMPOP=  -O3 -std=c11 -Wall
	#PROFILE= -O1 -p -g  
	COMPILER_FLAGS= -c $(COMPOP) $(PROFILE)
	BIND_FLAGS= -lm  $(PROFILE)
#   COMPILER_FLAGS= -c -s -O3 -funroll-loops
#   BIND_FLAGS= -s -O3  -funroll-loops
# -Zomf es para usar ilink.exe 5.0 que esta en \etc
#___________________________________________________________
#  flags for  Sun CC  (C++)  compiler
#   CC=CC 
#   COMPILER_FLAGS= -c -xO1 -w
#   BIND_FLAGS=
#___________________________________________________________
#  flags for  HP-UX ANSI C++ compiler (aCC, version A.03.50 an later) 
#   CC=aCC
#   COMPILER_FLAGS= -c +O1 -w -AA
#   BIND_FLAGS= -AA
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
    ./problemParam.h \
    ./prototypes.h

./advection.o: \
    ./advection.c \
    ./real.h \
    ./problemParam.h

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
    ./problemParam.h \
    ./real.h



./putfield.o: \
    ./putfield.c \
    ./real.h


#./getTakacsError.o: \
#    ./getTakacsError.c

./update.o: \
    ./update.c \
    ./real.h
    
./getDensity.o: \
    ./getDensity.c \
    ./real.h
    

#######################  do not change anything from here on ####################

    
%.o : %.c
	$(CC) $(COMPILER_FLAGS) $<

$(PROJ)$(EXT): $(objs)
	$(CC) $(objs)  -o  $@ $(BIND_FLAGS) $(COMPOP)
		
clean:
	@rm -rf $(objs)  *.o
	
