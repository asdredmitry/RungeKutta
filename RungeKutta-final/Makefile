help.o : help.c help.h
	gcc -o help.o -c -g help.c -lm
solve.o : solve.c solve.h
	gcc -o solve.o -c -g solve.c -lm
main.o : main.c
	gcc -o main.o -c -g  main.c -lm
main : main.o help.o solve.o
	gcc -g main.o help.o solve.o -lm -o main
clean : 
	rm -rf main.o help.o main solve.o
