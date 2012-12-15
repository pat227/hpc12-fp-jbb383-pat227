EXECUTABLES = householder wy BlockedQR BlockedQR2

all: $(EXECUTABLES)

householder: householder.c  matrices.o qrlib.o
	gcc -std=gnu99 -O3 -lrt -lm -Wall -Wextra  -o$@ $^
wy:	WYmain.c WY.o Utilities.o test.o MatrixMatrixMultiply.o MatrixTranspose.o
	gcc -std=gnu99 -Wall -Wextra -O3 -lrt -lm  -o$@ $^
BlockedQR:	BlockedQR_main.c BlockedQR.o Utilities.o WY.o MatrixTranspose.o
	gcc -fopenmp -std=gnu99 -Wall -Wextra -O3 -lrt -lm  -o$@ $^

BlockedQR2: BlockedQR_main2.c BlockedQR2.o Utilities.o WY.o MatrixTranspose.o BlockedQR.o
	gcc -fopenmp -std=gnu99 -Wall -Wextra -O3 -lrt -lm -o$@ $^

MatrixMatrixMultiply.o:	MatrixMatrixMultiply.c
	gcc -c -std=gnu99 -Wall -Wextra -Werror -g -o$@ $^

MatrixTranspose.o:	MatrixTranspose.c
	gcc -c -std=gnu99 -Wall -Wextra -Werror -g -o$@ $^
WY.o:	WY.c
	gcc -c -std=gnu99 -Wall -g -lm -o$@ $^
test.o: test.c
	gcc -c -std=gnu99 -Wall -lm  -g -o$@ $^
Utilities.o: Utilities.c
	gcc -c -lm -std=gnu99 -Wall -g -o$@ $^
BlockedQR.o: BlockedQR.c
	gcc -c -fopenmp -std=gnu99 -Wall -g -o$@ $^
BlockedQR2.o: BlockedQR2.c
	gcc -c -fopenmp -std=gnu99 -Wall -g -o$@ $^
matrices.o:	matrices.c
	gcc -c -std=gnu99 -lm -Wall -Wextra -Werror -o$@ $^
qrlib.o:	proj-LibCorrectness.c
	gcc -c -std=gnu99 -Wall -Wextra -Werror -o$@ $^
%.o: $.c %.h
	gcc $(DUBUGS_FLAGS) -c -std=gnu99 $<

clean:
	rm -f $(EXECUTABLES) *.o
