.PHONY=clean
CC=gcc
SOURCES=test.c br_tst.c fft_tst.c
OBJ_LIST=$(SOURCES:.c=.o)
CFLAGS= -W -Wall

PS6: $(OBJ_LIST)
	$(CC) $^ -o $@ -lm

debug: $(OBJ_LIST)
	$(CC) -g  $^ -o $@ 
.depend:
	gcc -MM $(SOURCES) > $@  

-include .depend

%.o:%.c
	$(CC) -o $@ -c $< $(CFLAGS) 
	
clean: 
	rm -f *.o ${OBJ_LIST}
