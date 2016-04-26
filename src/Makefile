# target : dep1 dep2 dep3
#	stmt1
#	stmt2
#	stmt3	
#
# $@ nom de la cible(target)
# $< nom de la premiere dependance
# $^  TOUTES les dependances

## juste des noms, pas des vrais fichiers
PROG=PBRandomReads
CFLAGS=-Wall -g

.PHONY:all clean

all: ../bin/${PROG}

../bin/${PROG} : ${PROG}.o utils.o
	mkdir -p ../bin
	gcc -o $@ $^

${PROG}.o : ${PROG}_draft.c utils.h
	gcc ${CFLAGS} -o $@ -c $<


utils.o : utils.c utils.h
	gcc ${CFLAGS} -o $@ -c $<

clean:
	rm -f ../bin/${PROG} *.o
