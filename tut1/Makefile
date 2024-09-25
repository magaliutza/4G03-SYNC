CFLAGS = -O2 --std=c++11
LFLAGS = 

Skeleton.T1: Skeleton.T1.o 
	g++ Skeleton.T1.o $(LFLAGS) -o Skeleton.T1

Skeleton.T1.o: Skeleton.T1.cc 
	g++ $(CFLAGS) -c Skeleton.T1.cc

clean:
	@\rm *.o Skeleton.T1
