include ../make.flags
include modules.flags

# Executables
all: res_evol.x  fil2fla.x jupdat2fla.x dat2fla.x dat2fla_clones.x orbit_diff.x vpladat2fla.x vastdat2fla.x detect_res.x #dmint_evol.x

LIBDEPS= ../lib/libprop.a ../lib/libsuit.a ../lib/libgauss.a ../lib/libmoid.a  

PROPNEO2FILES= pladat.o pla_dat.o planet_orbits.o read_pla.o angles_routines.o \
        rkg_local.o proutele.o prsingopt.o sortpre.o 

PRSINGEQFILES= prsingeq.o  chooserhs.o right_hand_side.o rad02pi.o \
      prdnod.o prextrdn.o prsalto.o rhsloop.o prstop.o prsector.o \
      falsectra.o choose_deg.o mutualrefcha.o kantorovich.o

DETRES= pladat.o pla_dat.o right_hand_side.o planet_orbits.o read_pla.o

#PRSINGEQFILES= prsingeq.o  chooserhs_chiara.o rhs2.o rad02pi.o \
      prdnod.o prextrdn.o prsalto.o rhsloop.o prstop.o prsector.o \
      falsectra.o choose_deg.o mutualrefcha.o kantorovich_chiara.o

RHSVAR=  pla_dat.o read_pla.o planet_orbits.o right_hand_side.o # rhs2.o

PLANETS=  pla_dat.o planet_orbits.o read_pla.o 

orbit_diff.x: orbit_diff.o $(LIBDEPS)
	$(FC) $(FFLAGS) orbit_diff.o $(LIBDEPS) -o orbit_diff.x

res_evol.x: res_evol.o  $(PROPNEO2FILES) $(PRSINGEQFILES) $(LIBDEPS) \
	../propag/quadp/libquamia.a 
	$(FC) $(FFLAGS) res_evol.o  $(PROPNEO2FILES) $(PRSINGEQFILES) \
	$(LIBDEPS) ../propag/quadp/libquamia.a  -o res_evol.x

dmint_evol.x: dmint_evol.o $(PROPNEO2FILES) $(PRSINGEQFILES) $(LIBDEPS) \
	dmintil_deriv.o
	$(FC) $(FFLAGS) dmint_evol.o $(PROPNEO2FILES) $(PRSINGEQFILES) \
	$(LIBDEPS) dmintil_deriv.o -o dmint_evol.x

#test_pla.x: test_pla.o $(PROPNEO2FILES) $(PRSINGEQFILES) $(LIBDEPS)
#	$(FC) $(FFLAGS) test_pla.o $(PROPNEO2FILES) $(PRSINGEQFILES) \
#	$(LIBDEPS) -o test_pla.x

fil2fla.x: fil2fla.o $(LIBDEPS)
	$(FC) $(FFLAGS) fil2fla.o $(LIBDEPS) -o fil2fla.x

dat2fla.x: dat2fla.o $(LIBDEPS)
	$(FC) $(FFLAGS) dat2fla.o $(LIBDEPS) -o dat2fla.x

jupdat2fla.x: jupdat2fla.o $(LIBDEPS)
	$(FC) $(FFLAGS) jupdat2fla.o $(LIBDEPS) -o jupdat2fla.x

vpladat2fla.x: vpladat2fla.o  $(LIBDEPS) 
	$(FC) $(FFLAGS) vpladat2fla.o $(LIBDEPS) -o vpladat2fla.x

vastdat2fla.x: vastdat2fla.o  $(LIBDEPS) 
	$(FC) $(FFLAGS) vastdat2fla.o $(LIBDEPS) -o vastdat2fla.x

dat2fla_clones.x: dat2fla_clones.o $(LIBDEPS)
	$(FC) $(FFLAGS) dat2fla_clones.o $(LIBDEPS) -o dat2fla_clones.x

detect_res.x: detect_res.o $(LIBDEPS) $(DETRES) ../propag/quadp/libquamia.a 
	$(FC) $(FFLAGS) detect_res.o $(LIBDEPS) $(DETRES) ../propag/quadp/libquamia.a -o detect_res.x


depend:
	../../bin/mkdep90 -I../include -L../suit -L../propag *.f90 *.h90>make.dep
	../../bin/mkmoddep -f ../.compiler ../suit ../propag > modules.flags

clean:
	-(rm -f *.x *.o *% *~ core *.mod)

distclean: clean
	../lib/deinst_src

neotar:	clean
	(cd ../ ; tar -cf neodys.tar neodys/)

# Dependencies
include make.dep

