all: tests

LIBEUTREEPATH=../euclidean_tree_lib

# ------------------------------------------------------------------- tests ---

tests: test_3d_strom #test_strom_z_clanku test_zlozitejsi_2d_strom

test_3d_strom: test_3d_strom.c
	gcc -c test_3d_strom.c -I$(LIBEUTREEPATH)
	gcc -o test_3d_strom test_3d_strom.o -lm -lgsl -lgslcblas -L$(LIBEUTREEPATH) -leutree

#test_strom_z_clanku:
#	echo "nada"
#
#test_zlozitejsi_2d_strom:
#	echo "nada"

# ------------------------------------------------------------------ servis ---

clean:
	rm -f *.o

deepclean: clean
	rm -f test_3d_strom
	rm -f test_strom_z_clanku
	rm -f test_zlozitejsi_2d_strom

