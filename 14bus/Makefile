# compiler option
CC = pgc++
NVCC = nvcc 

OPT = -acc -ta=nvidia -Minfo=accel 

#CUDA_PATH = /cm/extra/apps/CUDA/9.0.176_384.81/
#CUDA_PATH = /usr/local/cuda-9.0/

#INCLUDES := -I$(CUDA_PATH)/include
#LIBS := -L$(CUDA_PATH)/lib64/


SRCS = functions.cpp

#LIB = 

OBJ = LF14bus.o functions.o


EXE = lf14bus

run: all
	./lf14bus

all:
	$(CC) -acc -ta=nvidia -Minfo=accel LF14bus.cpp -c
	$(CC)  $(OPT) $(SRCS) $(LIB) -c 
	$(CC) $(OBJ) -o $(EXE) $(LIBS) $(OPT)  
clean: 
	rm -rf *.o ./lf14bus



