TARGET  = mpi_main

CC = mpicc
OPT =
CFLAGS = -Wall $(OPT)
LDFLAGS = #C:\Program Files (x86)\Microsoft SDKs\MPI\Include
OBJFILES = $(patsubst %.c,%.o,$(wildcard *.c))\
 $(patsubst %.c, ../serial_part/%.o, allocate_2D_array.c deallocate_2D_array.c)

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES)

%.o: %.c
	$(CC) $(CFLAGS)  -c -o $@ $<

.PHONY: clean
clean: 
	@echo "removing all object files"
	rm -rf *.o