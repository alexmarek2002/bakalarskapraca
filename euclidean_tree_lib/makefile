all: libeutree.a

# --------------------------------------------------------------------- lib ---

libeutree.a: et_console.o et_data.o et_geometry.o et_smith.o et_tikz.o
	ar -r libeutree.a et_console.o et_data.o et_geometry.o et_smith.o et_tikz.o

# -------------------------------------------------------------------- objs ---

et_console.o: et_console.c et_console.h et_data.h
	gcc -c et_console.c

et_data.o: et_data.c et_data.h
	gcc -c et_data.c

et_geometry.o: et_geometry.c et_geometry.h et_data.h
	gcc -c et_geometry.c

et_smith.o: et_smith.c et_smith.h et_data.h
	gcc -c et_smith.c

et_tikz.o: et_tikz.c et_tikz.h et_data.h
	gcc -c et_tikz.c

# ------------------------------------------------------------------ servis ---

clean:
	rm -f *.o

deepclean: clean
	rm -f libeutree.a

