EXECUTABLES = main 

all: $(EXECUTABLES)

main: main.c MatrixMatrixMultiply.o MatrixTranspose.o WY.o
	gcc -std=gnu99 -Wall -g -O3 $(DEBUG_FLAGS) -lrt -lm  -o$@ $^
MatrixMatrixMultiply.o:	MatrixMatrixMultiply.c
	gcc -c -std=gnu99 -Wall -g -o$@ $^
MatrixTranspose.o:	MatrixTranspose.c
	gcc -c -std=gnu99 -Wall -g -o$@ $^
WY.o:	WY.c
	gcc -c -std=gnu99 -Wall -g -lm -o$@ $^
%.o: $.c %.h
	gcc $(DUBUGS_FLAGS) -c -std=gnu99 $<

clean:
	rm -f $(EXECUTABLES) *.o
