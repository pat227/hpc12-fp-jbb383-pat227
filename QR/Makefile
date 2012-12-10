EXECUTABLES = main 

all: $(EXECUTABLES)

main: main.c MatrixMatrixMultiply.o MatrixTranspose.o WY.o test.o Utilities.o BlockedQR.o
	gcc -std=gnu99 -Wall -Wextra -g -O3 $(DEBUG_FLAGS) -lrt -lm  -o$@ $^
MatrixMatrixMultiply.o:	MatrixMatrixMultiply.c
	gcc -c -std=gnu99 -Wall -Wextra -Werror -g -o$@ $^
MatrixTranspose.o:	MatrixTranspose.c
	gcc -c -std=gnu99 -Wall -Wextra -Werror -g -o$@ $^
WY.o:	WY.c
	gcc -c -std=gnu99 -Wall -g -lm -o$@ $^
test.o: test.c
	gcc -c -std=gnu99 -Wall -lm  -g -o$@ $^
Utilities.o: Utilities.c
	gcc -c -std=gnu99 -Wall -g -o$@ $^

BlockedQR.o: BlockedQR.c
	gcc -c -std=gnu99 -Wall -g -o$@ $^

%.o: $.c %.h
	gcc $(DUBUGS_FLAGS) -c -std=gnu99 $<

clean:
	rm -f $(EXECUTABLES) *.o
