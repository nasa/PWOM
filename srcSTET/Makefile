include ../Makefile.def

SEARCH_EXTRA = -I${DATAREADINDICESDIR} -I${EMPIRICALIEDIR} -I${EMPIRICALUADIR} \
	-I${LIBDIR}

include ../Makefile.conf

include Makefile.DEPEND

include Makefile.RULES

MODULES = \
	ModGrid.o\
	ModSeMpi.o\
	ModMath.o\
	ModSePlot.o\
	ModBackground.o\
	ModSeState.o



METHODS = \

LIBOBJECTS = ${MODULES} ${METHODS}

EXEOBJECTS = main_stet.o ${LIBOBJECTS}

DEPEND:
	@${SCRIPTDIR}/depend.pl ${SEARCH} ${EXEOBJECTS} 


MY_LIB = ${LIBDIR}/libSTET.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo



${MY_LIB}: ${LIBOBJECTS}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${LIBOBJECTS}

EXE = ${BINDIR}/STET.exe

STET:	DEPEND
	@make ${EXE}
	@echo "STET.exe has been created"

${EXE}:	${EXEOBJECTS}
	@(if [ ${OS} = "Darwin" ];         then make STET_for_Darwin; \
	else                                    make STET_general; fi)

STET_general:
	${LINK.f90} -o ${EXE} ${EXEOBJECTS} \
		-L${LIBDIR}  -lEMPIRICALIE -lEMPIRICALUA -lINDICES -lSHARE -lTIMING \
		${Lflag1}

STET_for_Darwin:
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ${LIBDIR}/libEMPIRICALIE.a; \
		ar -x ${LIBDIR}/libEMPIRICALUA.a; \
		ar -x ${LIBDIR}/libINDICES.a; \
		ar -x ${LIBDIR}/libSHARE.a
	${LINK.f90} -o ${EXE} ${EXEOBJECTS} Tmp_/*.o ${LBLAS} ${Lflag1}
	rm -rf Tmp_

distclean: clean

