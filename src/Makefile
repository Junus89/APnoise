CC = gcc
CFLAGS = -g -Wall -Wno-unused-variable -std=c99
Lib = -lm
EX = APnoise

all : $(EX)
$(EX) : APnoise.o APnoisemain.o
	$(CC) $(CFLAGS) -o $@ $^ $(Lib)

BEMT.o : APnoise.c
	$(CC) $(CFLAGS) -o $@ -c $<

main.o : APnoisemain.c
	$(CC) $(CFLAGS) -o $@ -c $<

all:
	rm -rf ../bin
	mkdir ../bin && mv $(EX) ../bin

clean:
	rm -rf *.o rm -rf ../bin
