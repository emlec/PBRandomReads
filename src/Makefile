# Makefile structure

#target : dep1 dep2 dep3
#	stmt1
#	stmt2
#	stmt3

# $@ target name
# $< first dependency
# $^  all the dependencies


PROG = PBRandomReads_draft
CFLAGS = -Wall -g
CC = gcc

.PHONY:all clean mrproper

all: ../bin/${PROG}

../bin/${PROG} : ${PROG}.o utils.o
	mkdir -p ../bin
	${CC} ${CFLAGS} -o $@ $^

${PROG}.o: utils.h

%o:%c
	${CC} ${CFLAGS} -o $@ -c $<


clean:
	rm -f *.o

mrproper: clean
	rm -f ../bin/${PROG}
