CC=gcc -m64 
CFLAGS=-c -Wall -lm 

all: DMG

DMG: main.o angular.o slater.o file_io.o density.o strength.o
	$(CC) main.o angular.o slater.o file_io.o density.o strength.o -o DMG -lm -lgsl -lgslcblas


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

strength.o: strength.c
	$(CC) $(CFLAGS) strength.c


clean:
	rm -rf *.o DMG
