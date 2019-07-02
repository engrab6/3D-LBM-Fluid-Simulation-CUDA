

CFLAGS    =
FFLAGS    =
SOURCEC   = lbm.c
SOURCECXX = 
SOURCEF   =
SOURCEH   = lbm.h
DIRS      =
MANSEC    =
LOCDIR    =

-include Makefile.inc

include ${PETSC_DIR}/lib/petsc/conf/variables

test_lbm_OBJS = test_lbm.o
test_lbm_LIBS = -llbm

LIBBASE    = liblbm
LIBNAME    = ${LIBBASE}.${AR_LIB_SUFFIX}
CLEANFILES = ${LIBNAME} test_lbm

all: ${LIBNAME} test_lbm

include ${PETSC_DIR}/lib/petsc/conf/rules

${LIBNAME}: lib

test_lbm: ${test_lbm_OBJS} ${LIBNAME} chkopts
	-${CLINKER} -o $@ ${test_lbm_OBJS} ${PETSC_LIB} ${LIBNAME}
	-${DSYMUTIL} $@
	${RM} -f $(test_lbm_OBJS)