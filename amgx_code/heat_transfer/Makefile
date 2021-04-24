INC_D = -I/usr/local/cuda/include -I/home/divyaprakash/AMGX/base/include -I/home/divyaprakash/my_projects/fortarn_c_wrapper/amgx_code/heat_transfer_lib
LIB_D = -L/usr/local/cuda/lib64 -L/home/divyaprakash/AMGX/build -L/home/divyaprakash/my_projects/fortarn_c_wrapper/amgx_code/heat_transfer_lib
LIB_L = -lcusparse -lcublas -lcudart -ldl -lamgxsh -lht

solve_axb:
	gcc -O2 -std=c99 solve_axb.c -c $(INC_D)
	g++ -O2 solve_axb.o -o solve_axb $(LIB_D) $(LIB_L)
