CXX ?= g++
CFLAGS = -Wall -O3 -fPIC
SHVER = 2
OS = $(shell uname)

all: hhpssm

lib: svm.o
	if [ "$(OS)" = "Darwin" ]; then \
		SHARED_LIB_FLAG="-dynamiclib -W1,-install_name,libsvm.so.$(SHVER)"; \
	else \
		SHARED_LIB_FLAG="-shared -W1,-soname,libsvm.so.$(SHVER)"; \
	fi; \
	$(CXX) $${SHARED_LIB_FLAG} svm.o -o libsvm.so.$(SHVER)

hhpssm: src/hhpssm.cpp src/svm.o
	$(CXX) $(CFLAGS) -O3 src/hhpssm.cpp src/svm.o -o bin/hhpssm

clean:
	rm -f *~ bin/hhpssm src/svm.o src/libsvm.so.$(SHVER)
