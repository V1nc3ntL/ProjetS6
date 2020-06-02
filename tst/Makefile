.PHONY=clean
CC=gcc
SOURCES=test.c br_tst.c fft_tst.c fft_cmp_tst.c tst_func.c smp_tst.c cmp_tf.c
OBJ_LIST=$(SOURCES:.c=.o)
CFLAGS= -W -Wall

TEST: $(OBJ_LIST)
	$(CC) $^ -o $@   -lfftw3 -lm 

debug: $(OBJ_LIST)
	$(CC) -g  $^ -o $@ 
.depend:
	gcc -MM $(SOURCES) > $@  

-include .depend

%.o:%.c
	$(CC) -o $@ -c $< $(CFLAGS) 
	
clean: 
	rm -f *.o ${OBJ_LIST}
