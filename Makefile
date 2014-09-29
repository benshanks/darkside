APPS = mytest mytest2 mytest3
ROOTFLAGS = `root-config --cflags --glibs`
LIBFLAGS = -L$(ROOTSYS)/lib -lMathMore 



all: fitSpectrum
	g++ fitSpectrum.cxx -o fitSpectrum  $(LIBFLAGS) $(ROOTFLAGS)



clean:
	rm -f $(APPS) *.o *~ .depend

