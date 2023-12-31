TARGET = parallel_main

CC = gcc
OPT = -O2
CFLAGS = -Wall -fopenmp $(OPT) 
OBJFILES = $(patsubst %.c, %.o, $(wildcard *.c))\
 $(patsubst %.c, ../serial_part/%.o, allocate_2D_array.c deallocate_2D_array.c)

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean: 
	@echo "removing all objectfiles files"
	rm -rf *.o
print-%  : ; @echo $* = $($*)
