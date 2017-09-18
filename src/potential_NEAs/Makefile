include ../make.flags
include modules.flags

all:  planet_dist.x placros_NEAs.x mindist.x

#all:  possible_PNEAs.x get_PNEAs.x add_H.x propag_PNEAs.x \
	candidate_PNEAs.x planet_dist.x

# Executables


SECEVOL=secular_evolution.o dmintil_deriv.o right_hand_side.o 
libsec.a: $(SECEVOL)
	ar r ./libsec.a $(SECEVOL)
	$(RANLIB) ./libsec.a

PLANETS=planetorb_circ.o
libpla.a: $(PLANETS)
	ar r ./libpla.a $(PLANETS)
	$(RANLIB) ./libpla.a

#RHSVAR=  ../sec_evol/pla_dat.o ../sec_evol/read_pla.o \
	../sec_evol/planet_orbits.o ../sec_evol/right_hand_side.o 

LIBDEPS= ../lib/libprop.a ../lib/libmoid.a ../lib/libgauss.a ../lib/libsuit.a ./libsec.a ./libpert.a ../propag/quadp/libquamia.a ./libpla.a

#PROPNEO2= ../propneo2/pladat.o #../propneo2/clevpert.o 

PERTFUNC=perturbing_function.o
libpert.a: $(PERTFUNC)
	ar r ./libpert.a $(PERTFUNC)
	$(RANLIB) ./libpert.a

planet_dist.x: planet_dist.o spline1.o $(LIBDEPS)
	$(FC) $(FFLAGS) planet_dist.o spline1.o $(LIBDEPS) \
	-o planet_dist.x
mindist.x: mindist.o $(LIBDEPS)
	$(FC) $(FFLAGS) mindist.o $(LIBDEPS) -o mindist.x 

possible_PNEAs.x: possible_PNEAs.o $(LIBDEPS)
	$(FC) $(FFLAGS) possible_PNEAs.o $(LIBDEPS) -o possible_PNEAs.x

candidate_PNEAs.x: candidate_PNEAs.o  $(LIBDEPS) 
	$(FC) $(FFLAGS) candidate_PNEAs.o $(LIBDEPS) \
	 -o candidate_PNEAs.x

placros_NEAs.x: placros_NEAs.o  $(LIBDEPS) 
	$(FC) $(FFLAGS) placros_NEAs.o $(LIBDEPS) \
	 -o placros_NEAs.x

get_PNEAs.x: get_PNEAs.o
	$(FC) $(FFLAGS) get_PNEAs.o -o get_PNEAs.x

add_H.x: add_H.o
	$(FC) $(FFLAGS) add_H.o -o add_H.x

propag_PNEAs.x: propag_PNEAs.o $(LIBDEPS) spline1.o
	$(FC) $(FFLAGS) propag_PNEAs.o $(LIBDEPS) spline1.o -o propag_PNEAs.x

depend:
	../../bin/mkdep90 -I../include -L../suit -L../propag  *.f *.f90 *.h *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit ../propag  > modules.flags

clean:
	-(rm -f *.x *.o *% *~ *.mod core *.a a.out)


# Dependencies
include make.dep


