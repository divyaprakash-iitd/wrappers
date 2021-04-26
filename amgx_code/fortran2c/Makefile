INC_D = -I/home/divyaprakash/my_projects/fortarn_c_wrapper/amgx_code/heat_transfer_lib
LIB_D = -L/home/divyaprakash/my_projects/fortarn_c_wrapper/amgx_code/heat_transfer_lib
LIB_L = -lht

solve_axb:
	gfortran -c solveaxb.f90
	gfortran -o ftn_c.o solveaxb.o $(LIB_D) $(LIB_L)

