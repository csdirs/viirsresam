CXX=g++
LD=g++
CXXFLAGS=-g -O2 -Wall
LDFLAGS=-lhdf5 -lnetcdf -lz -lpthread -lm -lopencv_core
TARG=viirsresam
OFILES=\
	allocate_2d.o\
	readwrite.o\
	readwrite_ghrisst.o\
	create_viirs.o\
	resample.o\
	utils.o\
	main.o\

HFILES=\
	viirsresam.h\
	sort.h\

all: $(TARG)

$(TARG): $(OFILES)
	$(LD) -o $(TARG) $(OFILES) $(LDFLAGS)

%.o: %.cc $(HFILES)
	$(CXX) $(CXXFLAGS) -c $<

install: $(TARG)
	cp $(TARG) /usr/local/bin/

clean:
	rm -f $(OFILES) $(TARG)
