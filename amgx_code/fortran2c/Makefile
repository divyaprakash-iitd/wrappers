FORT= gfortran
LIB_D = -L/home/am/phd/amz208461/my_projects/wrappers/amgx_code/heat_transfer_lib
LIB_L = -lht


OBJECTS=ftn_c.o solveaxb.o

solveaxb: $(OBJECTS) 
	$(FORT) -o $@ $(OBJECTS) $(LIB_D) $(LIB_L)

ftn_c.o: ftn_c.f90 Makefile
	$(FORT) -c $<

solveaxb.o: solveaxb.f90 Makefile
	$(FORT) -c $<

clean:
	rm -f solveaxb *.mod *.bak *~ $(OBJECTS) 
