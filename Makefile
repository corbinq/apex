VERSION=0.3

CURL_STATIC_LIBS := $(shell curl-config --static-libs)
GSSAPI_STATIC_LIBS := $(shell krb5-config --libs gssapi)

CXXFLAGS= -static -std=c++11 -fopenmp -fno-math-errno -flto -O3 -pipe -s -fuse-ld=gold -I$(cdir)/src/eigen -I$(cdir)/src/spectra/include -I$(cdir)/src/boost -I$(cdir)/src/htslib/include -L$(cdir)/src/htslib/lib -I$(cdir)/src/BRENT -I$(cdir)/src/shrinkwrap/include
LDFLAGS= -static -lhts -fopenmp -flto -fuse-linker-plugin -O3 -pipe -Wno-unused-function -lgcov -lpthread -lz -lm -llzma -I$(cdir)/src/boost -I$(cdir)/src/htslib/include -L$(cdir)/src/htslib/lib -Wl,-rpath=$(cdir)/src/htslib/lib -I$(cdir)/src/BRENT -I$(cdir)/src/shrinkwrap/include -Wl,--no-as-needed
cppsrc = $(wildcard src/*.cpp src/BRENT/*.cpp)
csrc = $(wildcard src/*.c)

objs = $(cppsrc:.cpp=.o) $(csrc:.c=.o)

cdir = ${CURDIR}

bin/apex: $(objs)
	@mkdir -p $(@D)
	$(CXX) -o $@ $^ $(LDFLAGS)

.PHONY: clean

clean:
	rm -f *.o src/*.o src/BRENT/*.o


