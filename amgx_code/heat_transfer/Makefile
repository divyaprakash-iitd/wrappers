INC_D = -I/home/divyaprakash/my_projects/fortarn_c_wrapper/amgx_code/heat_transfer_lib
LIB_D = -L/home/divyaprakash/my_projects/fortarn_c_wrapper/amgx_code/heat_transfer_lib
LIB_L = -lht

solve_axb:
	gcc -O2 -std=c99 solve_axb.c -c $(INC_D)
	g++ -O2 solve_axb.o -o solve_axb $(LIB_D) $(LIB_L)
