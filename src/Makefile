# For Linux or any machines with gcc compiler
CC = gcc 
CFLAGS = -std=c99 -Wall -pedantic
BIN = ..\bin

all: Sec2_1 Sec2_3 Sec2_4

clean:
	del *.o *.exe $(BIN)\*.exe

OBJ = tiff.o allocate.o randlib.o qGGMRF.o solve.o

Sec2_1: Sec2_1.o $(OBJ)
	$(CC) $(CFLAGS) -o Sec2_1 Sec2_1.o $(OBJ) -lm

Sec2_3: Sec2_3.o $(OBJ) 
	$(CC) $(CFLAGS) -o Sec2_3 Sec2_3.o $(OBJ) -lm

Sec2_4: Sec2_4.o $(OBJ)
	$(CC) $(CFLAGS) -o Sec2_4 Sec2_4.o $(OBJ) -lm