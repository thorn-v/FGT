# Run
# ./bin/4gamete ./data/<file> > ./data/tst.out

# -O3 flag optimizes the code and makes a large difference in processing time
bin/4gamete: main.c
	gcc -O3 -o bin/4gamete main.c