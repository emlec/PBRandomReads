# Makefile structure

#target : dep1 dep2 dep3
#	stmt1
#	stmt2
#	stmt3

# $@ target name
# $< first dependency
# $^  all the dependencies

CC = gcc
CFLAGS = -Wall -g
LDFLAGS = -lz
EXEC = PBRandomReads_draft



.PHONY:all clean mrproper

all: ../bin/${EXEC}

../bin/${EXEC} : ${EXEC}.o utils.o fasta.o reads_distribution.o
	mkdir -p ../bin
	${CC} -o $@ $^ ${LDFLAGS}

${PROG}.o: utils.h fasta.h reads_distribution.h

%o:%c
	${CC} -o $@ -c $< ${CFLAGS}


clean:
	rm -f *.o *.gch

mrproper: clean
	rm -f ../bin/${EXEC}
