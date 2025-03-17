all: tests

# --------------------------------------------------------------------- lib ---

lib: et_console.o et_data.o et_geometry.o et_smith.o et_tikz.o

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

# ------------------------------------------------------------------- tests ---

tests: test_3d_strom test_strom_z_clanku test_zlozitejsi_2d_strom

test_3d_strom: lib test_3d_strom.c
	gcc -c test_3d_strom.c
	gcc -o test_3d_strom -lm -lgsl -lgslcblas test_3d_strom.o et_data.o et_console.o et_geometry.o et_smith.o et_tikz.o

test_strom_z_clanku:
	echo "nada"

test_zlozitejsi_2d_strom:
	echo "nada"

#snm: snm.o
#	gcc -o snm -lm -lgsl -lgslcblas snm.o
#
#snm.o: snm.c
#	gcc -c snm.c

# ------------------------------------------------------------------ servis ---

clean:
	rm -f *.o

deepclean: clean
	rm -f test_3d_strom
	rm -f test_strom_z_clanku
	rm -f test_zlozitejsi_2d_strom

