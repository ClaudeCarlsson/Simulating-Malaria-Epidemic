# Compiler (use mpicc which is the MPI C compiler wrapper)
CC = mpicc

# Compiler flags
CFLAGS = -Wall -Wextra -O3

# Math library
LIBS = -lm

# Name of the executable
EXE = simulate

# Source files
SRC = simulate.c prop.c

# Object files
OBJ = $(SRC:.c=.o)

# Default target
all: $(EXE)

# Link object files into the executable
$(EXE): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

# Compile source files into object files
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

# Clean up object files and the executable
.PHONY: clean
clean:
	rm -f $(OBJ) $(EXE)
