all: snm

snm: snm.o
	gcc -o snm -lm -lgsl -lgslcblas snm.o

snm.o: snm.c
	gcc -c snm.c

clean:
	rm -f *.o

deepclean: clean
	rm snm
