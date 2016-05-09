.PHONY: all clean

test all clean mrproper:
	(cd src && make $@)

