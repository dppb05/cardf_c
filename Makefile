CC=gcc
CFLAGS=-Wall -c
EXEDIR=.

all: $(EXEDIR)/cardf

$(EXEDIR)/cardf: stex.o matrix.o util.o cardf.o
	$(CC) stex.o matrix.o util.o cardf.o -o $@ -lm

cardf.o: cardf.c
	$(CC) $(CFLAGS) cardf.c

util.o: util.c
	$(CC) $(CFLAGS) util.c

matrix.o: matrix.c
	$(CC) $(CFLAGS) matrix.c

stex.o: stex.c
	$(CC) $(CFLAGS) stex.c

clean:
	rm -rf *o $(EXEDIR)/cardf
