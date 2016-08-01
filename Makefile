CC=gcc
CFLAGS=-Wall -O3 -c
EXEDIR=.

all: $(EXEDIR)/cardf

$(EXEDIR)/cardf: matrix.o util.o cardf.o
	$(CC) matrix.o util.o cardf.o -o $@ -lm

cardf.o: cardf.c
	$(CC) $(CFLAGS) cardf.c

util.o: util.c
	$(CC) $(CFLAGS) util.c

matrix.o: matrix.c
	$(CC) $(CFLAGS) matrix.c

clean:
	rm -rf *o $(EXEDIR)/cardf
