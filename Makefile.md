VERSION=0.3

CURL_STATIC_LIBS := $(shell curl-config --static-libs)
GSSAPI_STATIC_LIBS := $(shell krb5-config --libs gssapi)

CXXFLAGS= -std=c++11 -fopenmp -fno-math-errno -flto -O3 -pipe -s -fuse-ld=gold -I$(cdir)/eigen -I$(cdir)/spectra/include -I$(cdir)/boost -I$(cdir)/htslib/include -L$(cdir)/htslib/lib -I$(cdir)/BRENT -I$(cdir)/shrinkwrap/include
LDFLAGS= -fopenmp -flto -fuse-linker-plugin -O3 -pipe -Wno-unused-function -lgcov -lpthread -lz -lm -llzma -lhts -I$(cdir)/boost -I$(cdir)/htslib/include -L$(cdir)/htslib/lib -Wl,-rpath=$(cdir)/htslib/lib -I$(cdir)/BRENT -I$(cdir)/shrinkwrap/include -Wl,--as-needed
cppsrc = $(wildcard *.cpp BRENT/*.cpp)
csrc = $(wildcard *.c)

objs = $(cppsrc:.cpp=.o) $(csrc:.c=.o)

cdir = ${CURDIR}

bin/gqt: $(objs)
	@mkdir -p $(@D)
	$(CXX) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	rm -f *.o BRENT/*.o


