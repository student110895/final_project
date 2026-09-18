# Define compiler and strict compilation flags
CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
TARGET = symnmf

# Default target: compile the C program and link the math library (-lm)
$(TARGET): symnmf.c symnmf.h
	$(CC) $(CFLAGS) symnmf.c -lm -o $(TARGET)

# Clean target: remove the compiled executable from the directory
clean:
	rm -f $(TARGET)