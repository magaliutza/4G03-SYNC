CXX = c++
CFLAGS = -O2 -std=c++11

SQIsing: SQIsing.o MCSweeps.o
	$(CXX) $(CFLAGS) SQIsing.o MCSweeps.o -o SQIsing

SQIsing.o: SQIsing.cpp MCSweeps.h MCvar.h
	$(CXX) $(CFLAGS) -c SQIsing.cpp

MCSweeps.o: MCSweeps.cpp MCSweeps.h
	$(CXX) $(CFLAGS) -c MCSweeps.cpp
