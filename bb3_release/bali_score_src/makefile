install: bali_score 

clean:
	rm *.o

HEADERS = general.h clustalw.h

CC	= cc
CFLAGS  = -c -O -I$(EXPAT_INC)
LFLAGS	= -O -lm -L$(EXPAT_LIB) -lexpat
EXPAT_LIB	= /bach1/home/khanh2/lib
EXPAT_INC	= /bach1/home/khanh2/include


bali_score : readxml.o init.o util.o bali_score.o
	$(CC) -o $@ readxml.o init.o util.o bali_score.o $(LFLAGS)


.c.o :	$(HEADERS)
	$(CC) $(CFLAGS) $?

