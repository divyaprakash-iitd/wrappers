FORT= gfortran
LIB_D = -L/home/divyaprakash/my_projects/fortarn_c_wrapper/amgx_code/heat_transfer_lib -L/usr/lib/x86_64-linux-gnu 
LIB_L = -lht -lnvToolsExt



OBJECTS=ftn_c.o nvtx.o amgx_ht2d.o

amgx_ht2d: $(OBJECTS) 
	$(FORT) -o $@ $(OBJECTS) $(LIB_D) $(LIB_L)

ftn_c.o: ftn_c.f90 Makefile
	$(FORT) -c $<

nvtx.o: nvtx.f90 Makefile
	$(FORT) -c $<

amgx_ht2d.o: amgx_ht2d.f90 Makefile
	$(FORT) -c $<

clean:
	rm -f amgx_ht2d *.mod *.bak *~ $(OBJECTS) 
