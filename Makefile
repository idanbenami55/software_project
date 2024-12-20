CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LDFLAGS = -lm

symnmf: symnmf.o matrixutils.o symnmf.h
	$(CC) -o symnmf symnmf.o matrixutils.o $(LDFLAGS)

symnmf.o: symnmf.c matrixutils.h
	$(CC) -c symnmf.c $(CFLAGS)

matrixutils.o: matrixutils.c matrixutils.h
	$(CC) -c matrixutils.c $(CFLAGS)

clean:
	rm -f symnmf *.o 



