#make main

# -Wno-c++11-extensions —— turns off warnings about recent c++ features

CPP = clang++
CPPFLAGS = -c -Wno-c++11-extensions 

main: main.o laplace.o surface.o utilities.o
	$(CPP) main.o laplace.o surface.o utilities.o -o main.out
	rm *.o

main.o: main.cpp
	$(CPP) $(CPPFLAGS) main.cpp


laplace.o: laplace.cpp laplace.h
	$(CPP) $(CPPFLAGS) laplace.cpp


surface.o: surface.cpp surface.h
	$(CPP) $(CPPFLAGS) surface.cpp

utilities.o: utilities.cpp utilities.h
	$(CPP) $(CPPFLAGS) utilities.cpp


clean:
	rm *.o *~

#when name is different from makefile:
#make -f makefile main
