FORT= gfortran
LIB_D = -L/home/divyaprakash/my_projects/fortarn_c_wrapper/amgx_code/heat_transfer_lib
LIB_L = -lht


OBJECTS=ftn_c.o ht1dn.o

ht1dn: $(OBJECTS) 
	$(FORT) -o $@ $(OBJECTS) $(LIB_D) $(LIB_L)

ftn_c.o: ftn_c.f90 Makefile
	$(FORT) -c $<

ht1dn.o: ht1dn.f90 Makefile
	$(FORT) -c $<

clean:
	rm -f ht1dn *.mod *.bak *~ $(OBJECTS) 
