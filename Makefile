all: main

main: main.o normalize.o transform.o gauss_newton.o diff.o
	g++ -g main.o normalize.o transform.o gauss_newton.o diff.o -o main

main.o: main.cpp 
	g++ -g -c main.cpp -o main.o

normalize.o: normalize.cpp normalize.hpp
	g++ -g -c normalize.cpp -o normalize.o

transform.o: transform.cpp transform.hpp 
	g++ -g -c transform.cpp -o transform.o

gauss_newton.o: gauss_newton.cpp gauss_newton.hpp 
	g++ -g -c gauss_newton.cpp -o gauss_newton.o

diff.o: diff.cpp diff.hpp
	g++ -g -c diff.cpp -o diff.o

clean:
	rm *.o

