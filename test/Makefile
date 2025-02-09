CC=			gcc
CFLAGS=		-std=c99 -g -Wall -O3
CPPFLAGS=
INCLUDES=	-I..
OBJS=		msais32.o msais64.o
PROG=		msais-test
LIBS=		-lz

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

.SUFFIXES:.c .cpp .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

msais32.o:../msais.c ../msais.h
		$(CC) -c -o $@ $(CFLAGS) $<

msais64.o:../msais.c ../msais.h
		$(CC) -c -o $@ $(CFLAGS) -DMSAIS64 $<

msais-test:main.o $(OBJS)
		$(CC) -o $@ $< $(OBJS) $(LIBS)

clean:
		rm -fr *.o ../*.o a.out $(PROG) *~ *.a *.dSYM
