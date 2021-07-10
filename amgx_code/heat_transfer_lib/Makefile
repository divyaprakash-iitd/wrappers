# CUDA Toolkit location
INC_D = -I/usr/local/cuda/include -I/home/am/phd/amz208461/AMGX/base/include
LIB_D = -L/usr/local/cuda/lib64 -L/home/am/phd/amz208461/AMGX/build
LIB_L = -lcusparse -lcublas -lcudart -ldl -lamgxsh

heat_transfer:
	gcc -O2 -std=c99 -fPIC heat_transfer.c -c $(INC_D)
	g++ -O2 -shared heat_transfer.o $(LIB_D) $(LIB_L) -o libht.so
