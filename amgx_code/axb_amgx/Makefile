# CUDA Toolkit location
INC_D = -I/usr/local/cuda/include -I/home/am/phd/amz208461/AMGX/base/include
LIB_D = -L/usr/local/cuda/lib64 -L/home/am/phd/amz208461/AMGX/build
LIB_L = -lcusparse -lcublas -lcudart -ldl -lamgxsh

axb_amgx:
	gcc -O2 -std=c99 -fPIC axb_amgx.c -c $(INC_D)
	g++ -O2 -shared axb_amgx.o $(LIB_D) $(LIB_L) -o libaxbamgx.so
