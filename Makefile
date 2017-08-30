CC=g++
CFLAGS=-fopenmp -Wall -O3 -std=c++11 #-parallel -xHost -prof-gen -prof-dir=./  -g -fno-omit-frame-pointer -Wl,--no-as-needed -lprofiler -Wl,--as-needed

BASESRC=src/alphabet.cpp src/sequencefamily.cpp src/spacedworddb.cpp src/spacedfamily.cpp src/familyscore.cpp src/sequence.cpp src/patternset.cpp src/pattern.cpp
BASEHDR=src/alphabet.hpp src/sequencefamily.hpp src/spacedworddb.hpp src/spacedfamily.hpp src/familyscore.hpp src/sequence.hpp src/patternset.hpp src/spacedword.hpp src/pattern.hpp
BASEOBJ=$(BASESRC:.cpp=.o)

SWDBSRC=src/swdb.cpp src/dboptions.cpp src/createdatabase.cpp
SWDBHDR=src/dboptions.hpp src/createdatabase.hpp
SWDBOBJ=$(SWDBSRC:.cpp=.o)
SWDBPRG=swdb

SWDSSRC=src/swds.cpp src/dsoptions.cpp src/searchdatabase.cpp src/spacedhit.cpp
SWDSHDR=src/dsoptions.hpp src/searchdatabase.hpp src/spacedhit.hpp
SWDSOBJ=$(SWDSSRC:.cpp=.o)
SWDSPRG=swds

all: swdatabase swdatasearch

swdatabase: $(SWDBSRC) $(BASESRC) $(SWDBPRG)
swdatasearch: $(SWDSSRC) $(BASESRC) $(SWDSPRG)

$(SWDBPRG): $(SWDBOBJ) $(BASEOBJ)
	$(CC) $(CFLAGS) $(SWDBOBJ) $(BASEOBJ) -o $@
$(SWDSPRG): $(SWDSOBJ) $(BASEOBJ)
	$(CC) $(CFLAGS) $(SWDSOBJ) $(BASEOBJ) -o $@

.cpp.o: $(SWDBHDR) $(SWDSHDR) $(BASEHDR)
		$(CC) -c $(CFLAGS) $< -o $@

clean:
	find ./src/ -name "*.o" -delete
	find ./ -name $(SWDBPRG) -delete
	find ./ -name $(SWDSPRG) -delete

