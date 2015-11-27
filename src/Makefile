CC=g++

CFLAGS  = -O3 -m64 -std=c++0x -pthread 
LFLAGS  = -O3 -m64 -std=c++0x -lpthread

OMPFLAGS = -fopenmp


all: cut-ref processX VCF2FASTA-d VCF2FASTA-h VCF2VCFmin VCF2VDBV-d VCF2VDBV-h VDBV2VCFmin-d VDBV2VCFmin-h VDBV2FASTA tgc tgc_db

processX: processX.cpp
	$(CC) $(CFLAGS) processX.cpp -o $@

cut-ref: cut-ref.cpp
	$(CC) $(CFLAGS) cut-ref.cpp -o $@

VCF2FASTA-d: VCF2FASTA-d.cpp
	$(CC) $(CFLAGS) VCF2FASTA-d.cpp -o $@

VCF2FASTA-h: VCF2FASTA-h.cpp      
	$(CC) $(CFLAGS) VCF2FASTA-h.cpp -o $@  

VCF2VCFmin: VCF2VCFmin.cpp
	$(CC) $(CFLAGS) VCF2VCFmin.cpp -o $@  

VCF2VDBV-d: VCF2VDBV-d.cpp      
	$(CC) $(CFLAGS) $(OMPFLAGS) VCF2VDBV-d.cpp -o $@  

VCF2VDBV-h: VCF2VDBV-h.cpp        
	$(CC) $(CFLAGS) $(OMPFLAGS) VCF2VDBV-h.cpp -o $@ 

VDBV2VCFmin-d: VDBV2VCFmin-d.cpp        
	$(CC) $(CFLAGS) VDBV2VCFmin-d.cpp -o $@ 

VDBV2VCFmin-h: VDBV2VCFmin-h.cpp               
	$(CC) $(CFLAGS) VDBV2VCFmin-h.cpp -o $@ 

VDBV2FASTA: VDBV2FASTA.cpp
	$(CC) $(CFLAGS) VDBV2FASTA.cpp -o $@

.cpp.o:
	$(CC) $(CFLAGS) -c $<

tgc: tgc.o rangecod.o qsmodel.o
	$(CC) $(LFLAGS) -o $@ tgc.o rangecod.o qsmodel.o

tgc_db: tgc_db.o rangecod.o qsmodel.o
	$(CC) $(LFLAGS) -o $@ tgc_db.o rangecod.o qsmodel.o

clean:
	-rm cut-ref
	-rm processX
	-rm VCF2FASTA-d
	-rm VCF2FASTA-h
	-rm VCF2VCFmin
	-rm VCF2VDBV-d
	-rm VCF2VDBV-h
	-rm VDBV2VCFmin-d
	-rm VDBV2VCFmin-h
	-rm VDBV2FASTA
	-rm *.o
	-rm tgc
	-rm tgc_db