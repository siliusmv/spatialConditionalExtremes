GCC = gcc-12
spde_helpers = smat-operations.c smat-operations.h spde-precision.c spde-precision.h

.PHONY: test clean all

all: a.so b.so

a.so: a.c cgeneric.h
	$(GCC) -Wall -fpic -O -c -o a.o a.c -lm
	$(GCC) -shared -o a.so a.o -lm

b.so: b.c cgeneric.h $(spde_helpers)
	$(GCC) -Wall -fpic -O -c -o b.o b.c -lm
	$(GCC) -shared -o b.so b.o -lm

%.o: %.c %.h
	$(GCC) -Wall -fpic -O -c -o $@ $< -lm

clean:
	rm *.o
	rm *.so
