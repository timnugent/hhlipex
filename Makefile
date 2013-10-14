CPP = g++
CFLAGS = -Wall -Wextra -Werror -O3 -fPIC

all: hhpssm

src/svm.o: src/svm.cpp
	$(CPP) $(CFLAGS) -c src/svm.cpp -o src/svm.o

hhpssm: src/hhpssm.cpp src/svm.o
	$(CPP) $(CFLAGS) src/hhpssm.cpp src/svm.o -o bin/hhpssm

clean:
	rm bin/hhpssm src/*.o

test:
	./hhlipex.pl input/1gzmA.fa	
