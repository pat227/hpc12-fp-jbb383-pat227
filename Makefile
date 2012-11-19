EXECUTABLES = tests
#gcc -std=gnu99 -g -lrt -Wall -Wextra -Werror -o$@ $^
#gcc -std=gnu99 -O3 -march=native -mtune=native -ftree-vectorize -lrt -Wall -Wextra -Werror -o$@ $^
all: $(EXECUTABLES)

tests:	qr-test.c qrlib.o
	gcc -std=gnu99 -lcunit -Wall -Wextra -Werror -o$@ $^
qrlib.o:	proj-LibCorrectness.c
	gcc -c -std=gnu99 -Wall -Wextra -Werror -o$@ $^
clean:
	rm -f $(EXECUTABLES) *.o