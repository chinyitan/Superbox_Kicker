CC            = g77

LDFLAGS       = 

MAKEFILE      = Makefile

OBJS          = \
		convert.o \
		createna.o \
		disptens.o \
		dostop.o \
		fft.o \
		four1.o \
		geth.o \
		getpxx.o \
		getflags.o \
	        getrho.o \
		jacobi.o \
		makeext.o \
		makename.o \
		master.o \
		newcms.o \
		out_error.o \
		out_time.o \
		outdata.o \
		pusher1.o \
		pusher2.o \
		ppara.o \
		rdata.o \
		ritens.o \
		set_scales.o \
		sort.o \
		wdata.o \
		wsdata.o \
		ludcmp.o \
		lubksb.o \
	        super.o

PROGRAM       = super

SRCS_rec      = geth.recipe.f \
		getpxx.recipe.f 

SRCS_sup      = geth.f \
		getpxx.f 

SRCS          = \
		$(SRCS_sup) \
		convert.f \
		createna.f \
		disptens.f \
		dostop.f \
		fft.f \
		four1.f \
		getflags.f \
	        getrho.F \
		jacobi.f \
		makeext.f \
		makename.f \
		master.F \
		newcms.f \
		out_error.f \
		out_time.f \
		outdata.f \
		pusher1.f \
		pusher2.f \
		ppara.F \
		rdata.f \
		ritens.f \
		set_scales.f \
		sort.f \
		wdata.f \
		wsdata.f \
		ludcmp.f \
		lubksb.f \
	        super.F

norecipe:	
		cp -f geth.super.f geth.f 
		cp -f getpxx.super.f getpxx.f 
		cp -f fft.super.f fft.f 
recipe:	
		cp -f geth.recipe.f geth.f 
		cp -f getpxx.recipe.f getpxx.f 
		cp -f realft.f fft.f 

$(PROGRAM):     $(OBJS)  $(MAKEFILE) 
	        @echo "Linking $(PROGRAM) ..." 
		@$(FC) $(LDFLAGS)  $(OBJS)  -o $(PROGRAM).x.$(mesh)
	        @echo "done"

pgf90:  
	$(MAKE) "FC = pgf90" "FFLAGS = -O3 -D FORT90"  $(PROGRAM) 

pgf: 
	$(MAKE) "FC = pgf77" "FFLAGS = -O3"  $(PROGRAM)

g77: 
	$(MAKE) "FC = g77" "FFLAGS = -O6 -D GNUF77"   $(PROGRAM)

gfortran: 
	$(MAKE) "FC = gfortran" "FFLAGS = -O6 -D GNUF77 -Ofast -mcmodel=medium"   $(PROGRAM)

all:            
	$(PROGRAM)

.f.o:
convert.o: super.cb
createna.o: super.cb
disptens.o: super.cb 
dostop.o: super.cb 
eval.o: super.cb
fft.o: super.cb
four1.o: super.cb
getflags.o: super.cb
geth.o: super.cb 
geth.recipe.o: super.cb 
getrho.o: super.cb
getpxx.o: super.cb
getpxx.recipe.o: super.cb 
jacobi.o: super.cb 
makeext.o: super.cb
makename.o: super.cb
master.o: super.cb 
newcms.o: super.cb 
out_error.o: super.cb
out_time.o: super.cb
outdata.o: super.cb
pusher1.o: super.cb
pusher2.o: super.cb
ppara.o: super.cb
rdata.o: super.cb
ritens.o: super.cb 
set_scales.o: super.cb
sort.o: super.cb
wdata.o: super.cb 
wsdata.o:  super.cb
ludcmp.o: super.cb 
lubksb.o: super.cb 
super.o: super.cb


clean: 
	rm *.o 

#FFLAGS        = -O6 -m486 
#FFLAGS        = -g -pg -ansi -pedantic -Wall






