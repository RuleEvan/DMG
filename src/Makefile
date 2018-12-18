CC=gcc -m64
CFLAGS=-c -Wall -lm -ldl

all: SpeED-DMG

SpeED-DMG: main.o angular.o slater.o file_io.o density.o 
	$(CC) main.o angular.o slater.o file_io.o density.o -o SpeED-DMG -lm -ldl -lgsl -lgslcblas

main.o: main.c
	$(CC) $(CFLAGS) main.c

slater.o: slater.c
	$(CC) $(CFLAGS) slater.c

angular.o: angular.c
	$(CC) $(CFLAGS) angular.c

file_io.o: file_io.c
	$(CC) $(CFLAGS) file_io.c

density.o: density.c
	$(CC) $(CFLAGS) density.c

clean:
	rm -rf *.o SpeED-DMG
