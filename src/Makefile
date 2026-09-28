# conda 环境根目录
PREFIX     ?= $(CONDA_PREFIX)
INCLUDEDIR = $(PREFIX)/include
LIBDIR     = $(PREFIX)/lib

CPP      = $(CXX)                    
CPPFLAGS = -O3 -I$(PREFIX)/include   

LIBS = -Wl,-rpath,$(LIBDIR) -L$(LIBDIR) \
       -lboost_program_options -lhts -lz -llzma -ldeflate -lbz2 -lpthread -ldl

STDF       = -std=c++20
MYCFLAGS   = $(STDF) $(CPPFLAGS)

PROGS      = msavc_fasta vcf_merge
OBJMSAVC   = Arguments.o Fasta.o MultipleAlignmentFormat.o Mutation.o OutPut.o

all : $(PROGS)
	@echo done.

msavc_fasta: $(OBJMSAVC)
	$(CPP) msavc.cpp -o $@ $(OBJMSAVC) $(MYCFLAGS) $(LIBS)

vcf_merge: vcf_merge.cpp
	$(CPP) vcf_merge.cpp -o $@ $(MYCFLAGS) $(LIBS)

Arguments.o: Arguments.cpp Arguments.hpp Fasta.hpp MultipleAlignmentFormat.hpp
	$(CPP) $(MYCFLAGS) -c Arguments.cpp

Fasta.o: Fasta.cpp Fasta.hpp
	$(CPP) $(MYCFLAGS) -c Fasta.cpp

MultipleAlignmentFormat.o: MultipleAlignmentFormat.cpp MultipleAlignmentFormat.hpp Fasta.hpp
	$(CPP) $(MYCFLAGS) -c MultipleAlignmentFormat.cpp

Mutation.o: Mutation.cpp Mutation.hpp Fasta.hpp MultipleAlignmentFormat.hpp
	$(CPP) $(MYCFLAGS) -c Mutation.cpp

OutPut.o: OutPut.cpp OutPut.hpp Mutation.hpp Fasta.hpp MultipleAlignmentFormat.hpp
	$(CPP) $(MYCFLAGS) -c OutPut.cpp

msavc.o: msavc.cpp OutPut.hpp Arguments.hpp Fasta.hpp MultipleAlignmentFormat.hpp Mutation.hpp
	$(CPP) $(MYCFLAGS) -c msavc.cpp

clean:
	rm -f *.o *~ $(PROGS)