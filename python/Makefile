PYVER := $(shell script/config.py python_version)
OPTFLAGS := $(shell script/config.py opt_flags)
NUMPYINCLUDE := $(shell script/config.py numpy_include)
# readlink just normalizes the path:
LIB_PATH = $(shell readlink -m `pwd`/../lib)

CYTHON := cython
CXX := g++
LIBS := -lhermes2d-gcc-O-real -lglut -lJudy -lumfpack -lamd -lblas
CFLAGS := -Wfatal-errors -I/usr/include/python$(PYVER) $(OPTFLAGS) -I$(NUMPYINCLUDE)/numpy -I../include
LDFLAGS := -Wfatal-errors -Wl,--rpath=../lib:$(LIB_PATH) -L../lib/

HERMES := ../lib/libhermes2d-gcc-O-real.so

all: $(HERMES) hermes2d.so

$(HERMES):
	cd ../src; \
	make shared

examples: all c02.so c04.so c05.so c06.so c07.so c08.so

%.c: %.pyx
	$(CYTHON) $<

%.o : %.c
	$(CXX) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

%.so: %.o
	$(CXX) $(LDFLAGS) -shared $+ -o $@ $(LIBS)

clean:
	rm -f *.so *.pyd *.o *.pyc hermes2d.c

check:
	./test tests
