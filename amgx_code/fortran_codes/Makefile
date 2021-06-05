FORT= gfortran

OBJECTS=dataht.o ht2d.o

ht2d: $(OBJECTS) 
	$(FORT) -o $@ $(OBJECTS)

dataht.o: dataht.f90 Makefile
	$(FORT) -c $<

ht2d.o: ht2d.f90 Makefile
	$(FORT) -c $<

clean:
	rm -f ht1dn *.mod *.bak *~ $(OBJECTS) 
