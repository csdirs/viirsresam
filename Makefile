CXX=g++
LD=g++
CXXFLAGS=-g -O2 -Wall
LDFLAGS=-lhdf5 -lz -lpthread -lm -lopencv_core
TARG=viirsresam
OFILES=\
	allocate_2d.o\
	readwrite.o\
	resample_geo.o\
	main.o\

HFILES=viirsream.h

all: $(TARG)

$(TARG): $(OFILES)
	$(LD) -o $(TARG) $(OFILES) $(LDFLAGS)

%.o: %.cc $(HFILES)
	$(CXX) $(CXXFLAGS) -c $<

install: $(TARG)
	cp $(TARG) /usr/local/bin/

clean:
	rm -f $(OFILES) $(TARG)
