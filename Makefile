CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
TARGET = symnmf

$(TARGET): symnmf.c symnmf.h
	$(CC) $(CFLAGS) symnmf.c -lm -o $(TARGET)

clean:
	rm -f $(TARGET)