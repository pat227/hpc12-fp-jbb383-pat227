EXECUTABLES = householder wy BlockedQR BlockedQR2 BlockedQR_8 BlockedQR2_8 BlockedQR2_8_scaled BlockedQR_8_scaled

all: $(EXECUTABLES)

householder: householder.c matrices.o qrlib.o Utilities.o
	gcc -std=gnu99 -O3 -lrt -lm -Wall -Wextra  -o$@ $^
wy:	WYmain.c WY.o Utilities.o test.o MatrixMatrixMultiply.o MatrixTranspose.o
	gcc -std=gnu99 -Wall -Wextra -O3 -lrt -lm  -o$@ $^
BlockedQR:	BlockedQR_main.c BlockedQR.o Utilities.o WY.o MatrixTranspose.o MatrixMatrixMultiply.o test.o
	gcc -fopenmp -std=gnu99 -Wall -Wextra -O3 -lrt -lm  -o$@ $^
BlockedQR2: BlockedQR_main2.c BlockedQR2.o Utilities.o WY.o MatrixTranspose.o BlockedQR.o MatrixMatrixMultiply.o test.o
	gcc -fopenmp -std=gnu99 -Wall -Wextra -O3 -lrt -lm -o$@ $^
BlockedQR_8: 	BlockedQR_8_main.c BlockedQR8.o Utilities.o WY.o MatrixTranspose.o MatrixMatrixMultiply.o test.o
	gcc -fopenmp -std=gnu99 -Wall -Wextra -O3 -lrt -lm  -o$@ $^
BlockedQR2_8: BlockedQR_8_main2.c BlockedQR28.o Utilities.o WY.o MatrixTranspose.o BlockedQR.o MatrixMatrixMultiply.o test.o
	gcc -fopenmp -std=gnu99 -Wall -Wextra -O3 -lrt -lm -o$@ $^
BlockedQR_8_scaled: BlockedQR_8_scaled_main.c BlockedQR8.o Utilities.o WY.o MatrixTranspose.o MatrixMatrixMultiply.o test.o
	gcc -fopenmp -std=gnu99 -Wall -Wextra -O3 -lrt -lm  -o$@ $^
BlockedQR2_8_scaled: BlockedQR_8_scaled_main2.c BlockedQR28.o Utilities.o WY.o MatrixTranspose.o BlockedQR.o MatrixMatrixMultiply.o test.o
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
BlockedQR8.o: BlockedQR_8.c
	gcc -c -fopenmp -std=gnu99 -Wall -g -o$@ $^
BlockedQR28.o: BlockedQR2_8.c
	gcc -c -fopenmp -std=gnu99 -Wall -g -o$@ $^
matrices.o:	matrices.c
	gcc -c -std=gnu99 -lm -Wall -Wextra -Werror -o$@ $^
qrlib.o:	proj-LibCorrectness.c
	gcc -c -std=gnu99 -Wall -Wextra -Werror -o$@ $^
%.o: $.c %.h
	gcc $(DUBUGS_FLAGS) -c -std=gnu99 $<

clean:
	rm -f $(EXECUTABLES) *.o
