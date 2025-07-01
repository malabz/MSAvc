ENVNAME = package
PREFIX = ~/miniconda3/envs/$(ENVNAME)
INCLUDEDIR = $(PREFIX)/include
LIBDIR = $(PREFIX)/lib
DESTDIR = 

CPP = g++
#CPP = iCPP
CPPFLAGS = -O3 -static

LIBS = -lboost_program_options -lc -lstdc++ -ldl -lhts -lz -llzma -ldeflate -lbz2

STDF = -std=c++2a

MYCFLAGS =  $(STDF) $(CPPFLAGS)

PROGS = msavc_fasta vcf_merge

OBJMSAVC = Arguments.o Fasta.o MultipleAlignmentFormat.o Mutation.o OutPut.o

# HEADER = 

all : $(PROGS)
	@echo done.

msavc_fasta: $(OBJMSAVC)
	$(CPP) msavc.cpp -o $@ $(OBJMSAVC) $(MYCFLAGS) -L $(LIBDIR) $(LIBS)

vcf_merge: vcf_merge.cpp
	$(CPP) vcf_merge.cpp -o $@ $(MYCFLAGS)

Arguments.o: Arguments.cpp Arguments.hpp
	$(CPP) $(MYCFLAGS) -c Arguments.cpp -I $(INCLUDEDIR)

Fasta.o : Fasta.cpp Fasta.hpp
	$(CPP) $(MYCFLAGS) -c Fasta.cpp

MultipleAlignmentFormat.o : MultipleAlignmentFormat.cpp MultipleAlignmentFormat.hpp
	$(CPP) $(MYCFLAGS) -c MultipleAlignmentFormat.cpp

Mutation.o : Mutation.cpp Mutation.hpp
	$(CPP) $(MYCFLAGS) -c Mutation.cpp

OutPut.o : OutPut.cpp OutPut.hpp
	$(CPP) $(MYCFLAGS) -c OutPut.cpp -I $(INCLUDEDIR)

msavc.o: msavc.cpp OutPut.hpp Arguments.hpp
	$(CPP) $(MYCFLAGS) -c msavc.cpp

clean :
	rm -f *.o  *~ $(PROGS)
