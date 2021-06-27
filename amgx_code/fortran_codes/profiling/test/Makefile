FORT= gfortran
LIB_D = -L/usr/lib/x86_64-linux-gnu 
LIB_L = -lnvToolsExt


OBJECTS=nvtx.o check_nvtx.o

check_nvtx: $(OBJECTS) 
	$(FORT) -o $@ $(OBJECTS) $(LIB_D) $(LIB_L)

nvtx.o: nvtx.f90 Makefile
	$(FORT) -c $<

check_nvtx.o: check_nvtx.f90 Makefile
	$(FORT) -c $<

clean:
	rm -f check_nvtx *.mod *.bak *~ $(OBJECTS) 
