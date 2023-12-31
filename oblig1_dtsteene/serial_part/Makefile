TARGET  = serial_main

CC = gcc
OPT = -O2
CFLAGS = -Wall $(OPT)
OBJFILES = $(patsubst %.c, %.o, $(wildcard *.c))

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean: 
	@echo "removing all object files"
	rm -rf *.o