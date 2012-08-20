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
FPROG  = f1f209.f
CPROG  = MyTest.C 
# objects 
OBJ    = f1f209.o F1F209.o Radiator.o eInclusiveCrossSection.o MyTest.o 
# libraries 
LIBS   = -f2c -lgfortran
# directories 
SDIR   = ./src
IDIR   = ./include
ODIR   = ./obj

all: $(EXEC) 

$(EXEC): $(OBJ)
	mkdir -p bin
	mkdir -p obj
	$(CC) $(OFLAGS) bin/$(EXEC) $(OBJ) $(LIBS) 
	mv $(OBJ) $(ODIR) 

Radiator.o: $(SDIR)/Radiator.C $(IDIR)/Radiator.h
	$(CC) $(CFLAGS) $(SDIR)/Radiator.C

eInclusiveCrossSection.o: $(SDIR)/eInclusiveCrossSection.C $(IDIR)/eInclusiveCrossSection.h
	$(CC) $(CFLAGS) $(SDIR)/eInclusiveCrossSection.C

F1F209.o: $(SDIR)/F1F209.C $(IDIR)/F1F209.h
	$(CC) $(CFLAGS) $(SDIR)/F1F209.C

f1f209.o: $(SDIR)/$(FPROG)
	$(FF) $(CFLAGS) $(FFLAGS) $(SDIR)/$(FPROG) 

MyTest.o: $(CPROG)
	$(CC) $(CFLAGS) $(CPROG)

.PHONY: clean

clean: 
	rm $(ODIR)/*.o bin/$(EXEC)  
