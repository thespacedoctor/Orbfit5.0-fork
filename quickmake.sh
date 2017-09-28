# This script is used while developing to re-build only the essential components of Orbfit

cd /Users/Dave/Orbfit5.0/src/propag
gfortran -O3 -I../include -Wl,-rpath,/Users/Dave/anaconda/envs/orbfit5.0/lib -I../suit -c ephem_prop.f90 -o ephem_prop.o
ar r ../lib/libprop.a propag_state.o force_model.o rmodel.o  yark_pert.o  non_grav.o runge_kutta_gauss.o ra15_mod.o cobweb.o multiple_sol.o least_squares.o obs_correl.o rms_twobody.o count_opp.o count_opp2.o multi_store.o ephem_prop.o plocbd.o fitsubs.o pred_obs.o close_app.o tp_trace.o surf_trace.o semi_linear.o cla_store.o ret_analysistp.o offlov_checktp.o eval_risk.o obssto.o virtual_impactor.o fsteph_earth.o quadp/dqags.o quadp/dqagse.o quadp/dqagp.o quadp/dqagpe.o quadp/dqelg.o quadp/dqk21.o quadp/dqpsrt.o quadp/dqagsc.o quadp/dqagsec.o quadp/dqagpc.o quadp/dqagpec.o quadp/dqelgc.o quadp/dqk21c.o quadp/dqpsrtc.o  orbit_elements.o gaussdeg8.o arc_control.o laplace_poincare.o dyn_param.o  force_sat.o spher_harm.o perturbations.o force9d.o massmod.o des_routines5.o
cd /Users/Dave/Orbfit5.0/src/orbfit
make clean
make
