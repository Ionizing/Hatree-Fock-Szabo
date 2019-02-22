FORTRAN := gfortran
FLAGS   := -Wall -Waliasing -Wextra
LDFLAGS :=

build/original.out: build/original.o
	$(FORTRAN) -o $@ $< $(FLAGS) $(LDFLAGS)

build/original.o: src/original.f95
	$(FORTRAN) -c -o $@ $< $(FLAGS)
