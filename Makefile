APPS = fitSpectrum fitField fitFieldBeta
ROOTFLAGS = `root-config --cflags --glibs`
LIBFLAGS = -L$(ROOTSYS)/lib -lMathMore 


all: $(APPS)

fitSpectrum: fitSpectrum.cxx
	g++ fitSpectrum.cxx -o fitSpectrum  $(LIBFLAGS) $(ROOTFLAGS)

fitField: fitField.cxx
	g++ fitField.cxx -o fitField  $(LIBFLAGS) $(ROOTFLAGS)

fitFieldBeta: fitFieldBeta.cxx
	g++ fitFieldBeta.cxx -o fitFieldBeta  $(LIBFLAGS) $(ROOTFLAGS)

clean:
	rm -f $(APPS) *.o *~ .depend

