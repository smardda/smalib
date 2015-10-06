FC = ifort
#FFLAGS = -fdefault-real-8 -fdefault-double-8 -O
#FFLAGS = -fdefault-real-8 -O
FFLAGS = -r8 -save -O
#DBG =
#FFLAGS = -fdefault-real-8 -fdefault-double-8 -g -fbacktrace -fbounds-check -fexceptions
#FFLAGS = -fdefault-real-8 -g -fbacktrace -fbounds-check -fexceptions
FFLAGS = -r8 -save -debug -debug-parameters all -trace -fpe0 -fp-stack-check -check
DBG = g

src =  bsplvn.f  bvalue.f  interu.f  interv.f  tb06a.f  tc01a.f  tg04a.f \
	isnumb.f randlx.f spleval.f msrect.f lowor.f \
	isort.f \
	cfftb1.f cffti1.f passb3.f passb5.f passf2.f passf4.f passf.f zfftf.f \
	cfftf1.f passb2.f passb4.f passb.f passf3.f passf5.f zfftb.f zffti.f


obj = $(src:.f=.o)

lib$(DBG).a: $(obj)
	ar ru lib$(DBG).a $(obj)
	ranlib lib$(DBG).a

tar: $(src) README makefile
	tar cvf fort.tar $(src) README makefile

clean: $(obj)
	rm -f $(obj)
