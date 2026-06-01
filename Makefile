CC = g++
CCFLAGS = -g

all: main

main: main.o normalize.o transform.o gauss_newton.o diff.o integration.o helper_functions.o chol.o
	$(CC) $(CCFLAGS) main.o normalize.o transform.o gauss_newton.o diff.o integration.o helper_functions.o chol.o -o main

main.o: main.cpp 
	$(CC) $(CCFLAGS) -c main.cpp -o main.o

normalize.o: normalize.cpp normalize.hpp
	$(CC) $(CCFLAGS) -c normalize.cpp -o normalize.o

transform.o: transform.cpp transform.hpp 
	$(CC) $(CCFLAGS) -c transform.cpp -o transform.o

gauss_newton.o: gauss_newton.cpp gauss_newton.hpp 
	$(CC) $(CCFLAGS) -c gauss_newton.cpp -o gauss_newton.o

diff.o: diff.cpp diff.hpp
	$(CC) $(CCFLAGS) -c diff.cpp -o diff.o

integration.o: integration.cpp integration.hpp
	$(CC) $(CCFLAGS) -c integration.cpp -o integration.o

helper_functions.o: helper_functions.cpp helper_functions.hpp
	$(CC) $(CCFLAGS) -c helper_functions.cpp -o helper_functions.o

chol.o: chol.cpp chol.hpp
	$(CC) $(CCFLAGS) -c chol.cpp -o chol.o

clean:
	rm -f *.o main
