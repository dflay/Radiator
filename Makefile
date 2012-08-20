# makefile for compiling and linking C++ and fortran 

# compiler 
FF     = gfortran
CC     = g++
# flags 
CFLAGS = -c
FFLAGS = -fno-underscoring
OFLAGS = -o 
# executable  
EXEC   = Test
# programs
CPROG  = MyTest.C 
# objects 
#OBJ    = f1f209.o F1F209.o Radiator.o eInclusiveCrossSection.o 

# object files to build which have corresponding header and c file
CPPOBJECTS := obj/Radiator.o \
              obj/eInclusiveCrossSection.o \
              obj/F1F209.o
# object files to build from fortran files
F77OBJECTS := obj/f1f209.o

OBJ := $(CPPOBJECTS) $(F77OBJECTS)

# libraries 
LIBS   = -lgfortran
# directories 
SDIR   = src
IDIR   = include
ODIR   = obj

all: $(EXEC) 

$(EXEC): $(CPPOBJECTS) $(F77OBJECTS) 
	mkdir -p bin
	mkdir -p obj
	$(CC) -o bin/$(EXEC) $(CPROG) $(OBJ) $(LIBS) 

$(CPPOBJECTS) : $(ODIR)/%.o : $(IDIR)/%.h
	$(CC) -o $@ $(CFLAGS) $(@:$(ODIR)/%.o=$(SDIR)/%.C)

$(F77OBJECTS) : $(ODIR)/%.o :
	$(FF) -o $@ $(CFLAGS) $(FFLAGS) $(@:$(ODIR)/%.o=$(SDIR)/%.f)

#Radiator.o: $(SDIR)/Radiator.C $(IDIR)/Radiator.h
#	$(CC) $(CFLAGS) $(SDIR)/Radiator.C
#
#eInclusiveCrossSection.o: $(SDIR)/eInclusiveCrossSection.C $(IDIR)/eInclusiveCrossSection.h
#	$(CC) $(CFLAGS) $(SDIR)/eInclusiveCrossSection.C
#
#F1F209.o: $(SDIR)/F1F209.C $(IDIR)/F1F209.h
#	$(CC) $(CFLAGS) $(SDIR)/F1F209.C
#
#f1f209.o: $(SDIR)/$(FPROG)
#	$(FF) $(CFLAGS) $(FFLAGS) $(SDIR)/$(FPROG) 
#
#MyTest.o: $(CPROG)
#	$(CC) $(CFLAGS) $(CPROG)

.PHONY: clean

clean: 
	rm $(ODIR)/*.o bin/$(EXEC) 


