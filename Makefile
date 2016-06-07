RM     = rm -f
SHELL  = /bin/sh

MYOS := $(subst -,,$(shell uname))

ifeq ($(MYOS),HPUX)
  IMSLIBS = /usr/site2/imsl/libimsl.a
#  FFLAGS=  -O -C  +es +Obb1000 -w -K -R8 +FPVZOU
#  FFLAGS=  -g -C  +es +Obb1000 -K -R8 +FPZOU
  FFLAGS=  -g -C  +es +Obb1000 -K +FPZOU
  LDFLAGS=-Wl,-a archive
  OTHERLIBS = -Wl,-L$(CODA)/HP_UX/lib \
	-lana -lmsg -lcoda -Wl,-L$(CERN_ROOT)/lib -lgenlib -lpacklib \
        -lkernlib
endif

ifeq ($(MYOS),ULTRIX)
  FFLAGS=-check_bounds
  LDFLAGS=
  OTHERLIBS = -L$(CODA)/ULTRIX/lib \
	-lana -lmsg -lcoda -L$(CERN_ROOT)/lib -lpacklib
endif

ifeq ($(MYOS),OSF1)
  FFLAGS=  -O -C -extend_source -r8 -fpe -warn unused
#  FFLAGS=  -O -C -extend_source -fpe
endif


#ifeq ($(MYOS),Linux)  
#   F77=g77
#   FFLAGS=-g -ffixed-line-length-132
#   FFLAGS=-g -I$(Csoft)/INCLUDE -ffixed-line-length-132
#   DISPFLAGS=$(FFLAGS)
#   OTHERLIBS += -lc -lm -lnsl
#   OURLIBS := $(OURGENLIBS) $(LIBROOT)/libport.a
#endif

ifeq ($(MYOS),Linux)
# Uncomment these two lines for redhat 7.2
#  ABSOFT=/apps/absoft/PRO/usr/absoft
#  CERN_ROOT=/apps/cernlib/i386_redhat72/2001
# Uncomment these two lines for Enterprise Linux 3.
  ABSOFT=/apps/absoft/absoft-8.2/opt/absoft/
  CERN_ROOT = /apps/cernlib/i386_rhel3/2003
  FABSFLAGS=-O -V -W -f -s -N1 -B108 -B100 -N90 -N22 -N2 -N113
  INCLUDES=-I.
  EXTRAFLAGS=-DABSOFTFORTRAN
  FFLAGS= $(INCLUDES) $(FABSFLAGS) $(EXTRAFLAGS)
  FFLAG1=$(FFLAGS) -c
  OTHERLIBS = -L$(LIBROOT) -lctp \
        -L$(CERN_ROOT)/lib $(CERNLIBS) -lV77 -lU77 -lg2c -lc -lm \
	-lnsl -lcrypt
# Uncomment the last line above (-lnsl -lcrypt) if you are using cernlib 2001.
  FC  := $(ABSOFT)/bin/f77
  F77 :=$(ABSOFT)/bin/f77
endif

radcalc_objs = dilog.o f_to_sig.o fy.o model_fix.o  bdisnew4he3.o w1w2.o\
              nform.o quadmo.o radiate_calc.o sigmodel_calc.o dis_smear.o target_info.o\
              xsechd.o y_calc.o evalspline.o cubgcv.o interv.o getkine.o rlt.o r1990.o \
              radcalc.o eval.o f2glob.o sig_bar_df.o f1f2in06.o atailfl.o

xem_radcalc_objs = dilog.o f_to_sig.o fy.o model_fix.o bdisnew4he3.o\
              nform.o quadmo.o radiate_calc.o sigmodel_calc.o dis_smear.o\
              w1w2.o target_info.o\
              xsechd.o y_calc.o evalspline.o cubgcv.o interv.o getkine.o rlt.o r1990.o\
              xem_radcalc.o eval.o f2glob.o sig_bar_df.o f1f2in06.o atailfl.o

xem_radcalc_ext_objs = dilog.o f_to_sig.o fy.o model_fix.o bdisnew4he3.o w1w2.o\
              nform.o quadmo.o radiate_calc_ext.o sigmodel_calc.o dis_smear.o target_info_ext.o\
              xsechd.o y_calc.o evalspline.o cubgcv.o interv.o getkine.o rlt.o r1990.o\
              xem_radcalc_ext.o eval.o f2glob.o sig_bar_df.o  f1f2in06.o atailfl.o

none: radcalc xem_radcalc xem_radcalc_external

all: radcalc xem_radcalc xem_radcalc_external

radcalc: $(radcalc_objs) Makefile
	 $(F77) -o $@ $(FFLAGS) $(radcalc_objs) $(OTHERLIBS) $(IMSLIBS)

xem_radcalc: $(xem_radcalc_objs) Makefile
	 $(F77) -o $@ $(FFLAGS) $(xem_radcalc_objs) $(OTHERLIBS) $(IMSLIBS)

xem_radcalc_external: $(xem_radcalc_ext_objs) Makefile
	 $(F77) -o $@ $(FFLAGS) $(xem_radcalc_ext_objs) $(OTHERLIBS) $(IMSLIBS)

clean: 
	 rm -f *.o radcalc xem_radcalc xem_radcalc_external
