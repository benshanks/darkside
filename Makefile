APPS = fitSpectrum fitField fitFieldBeta
ROOTFLAGS = `root-config --cflags --glibs`
LIBFLAGS = -L$(ROOTSYS)/lib -lMathMore 


all: $(APPS)

%: %.cxx
	g++ $< -o $@  $(LIBFLAGS) $(ROOTFLAGS)

clean:
	rm -f $(APPS) *.o *~ .depend

