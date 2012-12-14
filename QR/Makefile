EXECUTABLES = householder wy

all: $(EXECUTABLES)

#main: main.c MatrixMatrixMultiply.o MatrixTranspose.o WY.o test.o Utilities.o BlockedQR.o MatrixVectorMultiply.o
#	gcc -std=gnu99 -Wall -Wextra -g -O3 $(DEBUG_FLAGS) -lrt -lm  -o$@ $^
householder: householder.c  matrices.o qrlib.o
	gcc -std=gnu99 -lrt -lm -Wall -Wextra -Werror -o$@ $^
wy:	WYmain.c WY.o Utilities.o test.o MatrixMatrixMultiply.o
	gcc -std=gnu99 -Wall -Wextra -O3 -lrt -lm  -o$@ $^
MatrixMatrixMultiply.o:	MatrixMatrixMultiply.c
	gcc -c -std=gnu99 -Wall -Wextra -Werror -g -o$@ $^
#MatrixVectorMultiply.o:	MatrixVector.c
#	gcc -c -std=gnu99 -Wall -Wextra -Werror -g -o$@ $^
#MatrixTranspose.o:	MatrixTranspose.c
#	gcc -c -std=gnu99 -Wall -Wextra -Werror -g -o$@ $^
WY.o:	WY.c
	gcc -c -std=gnu99 -Wall -g -lm -o$@ $^
test.o: test.c
	gcc -c -std=gnu99 -Wall -lm  -g -o$@ $^
Utilities.o: Utilities.c
	gcc -c -lm -std=gnu99 -Wall -g -o$@ $^
BlockedQR.o: BlockedQR.c
	gcc -c -std=gnu99 -Wall -g -o$@ $^
matrices.o:	matrices.c
	gcc -c -std=gnu99 -lm -Wall -Wextra -Werror -o$@ $^
qrlib.o:	proj-LibCorrectness.c
	gcc -c -std=gnu99 -Wall -Wextra -Werror -o$@ $^
%.o: $.c %.h
	gcc $(DUBUGS_FLAGS) -c -std=gnu99 $<

clean:
	rm -f $(EXECUTABLES) *.o
