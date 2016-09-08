include ../../Makefile.def

SEARCH_EXTRA = -I../ -I${DATAREADINDICESDIR} -I${EMPIRICALIEDIR} -I${EMPIRICALUADIR} \
	-I${LIBDIR}

include ../../Makefile.conf

MPIRUN = mpirun -np 2

help:
	@echo "Tests for Particle module: unit tests and combinations"
	@echo
	@echo "   test_axes             - test CON_axes"
	@echo "   test_coord_transform  - test ModTimeConvert"
	@echo "   test                  - test all the above"
	@echo
	@echo "To test everything with easy to read output try"
	@echo
	@echo "   make -s test"
	@echo
	@echo "The parallel tests can be executed like this"
	@echo
	@echo "   make MPIRUN='mpirun -np 4' test_line_extract"
	@echo
	@echo "The default value for MPIRUN is mpirun -np 2."

test:   unit_test_sample
	@echo ---------------------------------------------------------
	@ls -l *.diff

distclean: clean
	rm -rf *.tmp *.diff test*.out test_initial_state.ref \
		test_process_var_name.log

#######################################################################
OBJ_SAMPLE = \
	../ModParticle.o 


unit_test_sample_exe:
	@(cd .. ;make ModParticle.o)
	@make unit_test_sample.o 
	${LINK.f90} -o unit_test_sample.exe ${OBJ_SAMPLE} \
		unit_test_sample.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_sample: 
	@echo unit_test_sample starting
	@make unit_test_sample_exe 
	@echo ./unit_test_sample_exe
	-@(./unit_test_sample_exe )

#######################################################################

OBJ_PUSHER = \
	../ModParticle.o 


unit_test_pusher_exe:
	@(cd .. ;make ModParticle.o)
	@make unit_test_pusher.o 
	${LINK.f90} -o unit_test_pusher.exe ${OBJ_PUSHER} \
		unit_test_pusher.o -L${LIBDIR}  -lEMPIRICALIE \
		-lEMPIRICALUA -lINDICES -lSHARE -lTIMING ${Lflag1}

unit_test_pusher: 
	@echo unit_test_pusher starting
	@make unit_test_pusher_exe 
	@echo ./unit_test_pusher_exe
	-@(./unit_test_pusher_exe )

#######################################################################

