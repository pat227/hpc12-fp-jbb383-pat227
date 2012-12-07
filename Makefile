EXECUTABLES = tests matrices_tests householder
#gcc -std=gnu99 -g -lrt -Wall -Wextra -Werror -o$@ $^
#gcc -std=gnu99 -O3 -march=native -mtune=native -ftree-vectorize -lrt -Wall -Wextra -Werror -o$@ $^
all: $(EXECUTABLES)

tests:	qr-test.c qrlib.o
	gcc -std=gnu99 -lcunit -Wall -Wextra -Werror -o$@ $^
qrlib.o:	proj-LibCorrectness.c
	gcc -c -std=gnu99 -Wall -Wextra -Werror -o$@ $^
matrices_tests:	matrices_tests.c matrices.o
	gcc -std=gnu99 -lcunit -lm -Wall -Wextra -Werror -o$@ $^
matrices.o:	matrices.c
	gcc -c -std=gnu99 -lm -Wall -Wextra -Werror -o$@ $^
householder: householder.c  matrices.o qrlib.o
	gcc -std=gnu99 -lm -Wall -Wextra -Werror -o$@ $^
clean:
	rm -f $(EXECUTABLES) *.o