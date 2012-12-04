EXECUTABLES = main 

all: $(EXECUTABLES)

main: main.c MatrixMatrixMultiply.c MatrixTranspose.c  #test.c
	gcc -std=gnu99 -Wall -g -O3 $(DEBUG_FLAGS) -lrt -o$@ $^

%.o: $.c %.h
	gcc $(DUBUGS_FLAGS) -c -std=gnu99 $<

clean:
	rm -f $(EXECUTABLES) *.o
