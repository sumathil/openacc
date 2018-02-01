
TARGET = lf14bus
CC = icc
CFLAGS = -std=c++11 -fopenmp -xCORE-AVX2 -axMIC-AVX512 -O3  

DEP = LF_include.h LF_functions.h
OBJ = LF14bus.o functions.o

%.o: %.cpp $(DEP)
	$(CC) -o $(CFLAGS)
lf14bus: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

all: $(TARGET)

run: all
	- ./lf14bus
clean:
	rm -rf *.o lf14bus	


