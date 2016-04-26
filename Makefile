#
# Syntaxe
#
# target : dep1 dep2 dep3
#	stmt1
#	stmt2
#	stmt3		 
#
.PHONY: all clean

## all OU clean
all clean:
	(cd src && make $@)

