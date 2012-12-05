EXECUTABLES = main 

all: $(EXECUTABLES)

main: main.c MatrixMatrixMultiply.c MatrixTranspose.c WY.c
	gcc -std=gnu99 -Wall -g -O3 $(DEBUG_FLAGS) -lrt -lm  -o$@ $^

%.o: $.c %.h
	gcc $(DUBUGS_FLAGS) -c -std=gnu99 $<

clean:
	rm -f $(EXECUTABLES) *.o
