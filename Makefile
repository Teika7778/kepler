all: main

main: main.o normalize.o transform.o
	g++ -g main.o normalize.o transform.o -o main

main.o: main.cpp 
	g++ -g -c main.cpp -o main.o

normalize.o: normalize.cpp normalize.hpp
	g++ -g -c normalize.cpp -o normalize.o

transform.o: transform.cpp transform.hpp 
	g++ -g -c transform.cpp -o transform.o

clean:
	rm *.o

