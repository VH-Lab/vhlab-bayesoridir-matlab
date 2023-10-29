

Notes:

To compile hdf5 on Apple Silicon, use -L/opt/homebrew/lib in addition to -lhdf5

For example:

gcc arrayhelp.c arrayhelp_test.c filehelp.c -O3 -L/opt/homebrew/lib -lhdf5 -o ./test
