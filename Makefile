FORTRAN := gfortran
FLAGS   := -Wall -Waliasing -Wextra
LDFLAGS :=

run: build/original.out
	build/original.out

build/original.out: build/original.o
	$(FORTRAN) -o $@ $< $(FLAGS) $(LDFLAGS)

build/original.o: src/original.f95
	mkdir -p build
	$(FORTRAN) -c -o $@ $< $(FLAGS)
