CPP=g++
F90=gfortran
CFLAGS  = -std=c++0x
FFLAGS  = -ffixed-line-length-none
LDFLAGS = -lgfortran

SOURCES_C = gb5.cpp func.cpp vector_class.cpp
SOURCES_F = defs_basis.F90 readin.F90 42_parser/instrng.F90 42_parser/incomprs.F90 42_parser/inreplsp.F90  32_util/inupper.F90 \
42_parser/intagm.F90 42_parser/inarray.F90 42_parser/inread.F90 32_util/isfile.F90 

OBJECTS_C=$(SOURCES_C:.cpp=.o)
OBJECTS_F=$(SOURCES_F:.F90=.o)
EXECUTABLE=hello

all: $(SOURCES_C) $(SOURCES_F) $(EXECUTABLE)
	mv gb5 ..

$(EXECUTABLE): $(OBJECTS_C) $(OBJECTS_F) 
	$(CC) $(LDFLAGS) $(OBJECTS_C) $(OBJECTS_F) -o $@

.cpp.o:
	$(CC)  $(CFLAGS) $< -o $@

.F90.o:
	$(F90) $(FFLAGS) $< -o $@

clean:
	rm *o
