CC = gcc
CFLAGS = -std=c99 -O2
UNIX = -lm

all: main

main: main.o dst.o debug_utils.o
	${CC} -g main.o dst.o debug_utils.o -o ex4 ${UNIX}
		
main.o: main.c
	${CC} ${CFLAGS} -c main.c ${UNIX}

dst.o: dst.c dst.h
	${CC} ${CFLAGS} -c dst.c ${UNIX}

debug_utils.o: debug_utils.c debug_utils.h
	${CC} ${CFLAGS} -c debug_utils.c ${UNIX}

clean:
	rm -rf bf *.o

.PHONY: clean
