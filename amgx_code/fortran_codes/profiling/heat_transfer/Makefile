FORT= f95
LIB_D = -L/home/apps/centos7/nvidia-hpc-sdk/20.11/cuda11.0/Linux_x86_64/20.11/cuda/lib64 -L/home/am/phd/amz208461/my_projects/wrappers/amgx_code/heat_transfer_lib
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
