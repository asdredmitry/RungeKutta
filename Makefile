help.o : help.c help.h
	gcc -o help.o -c help.c -lm
solve.o : solve.c solve.h
	gcc -o solve.o -c solve.c -lm
main.o : main.c
	gcc -o main.o -c main.c -lm
main : main.o help.o solve.o
	gcc main.o help.o solve.o -lm -o main
clean : 
	rm -rf main.o help.o main solve.o