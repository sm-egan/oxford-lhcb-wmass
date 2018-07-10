ROOT=`root-config --cflags` `root-config --libs --glibs` 

#LHAPDF_DIR=/home/pili/local

#LHAPDF=-I${#LHAPDF_DIR}/include -L${#LHAPDF_DIR}/lib/ -lLHAPDF
#CFLAGS = -O3  -Wall -std=c++14 -pedantic -Wshadow
CFLAGS = -O3  -Wall -std=c++1y -pedantic -Wshadow
#CC=g++
CC=$(CXX)

all: bin/create_templates.exe

bin/%.exe: src/%.cxx
	@mkdir -p tmp/
	@mkdir -p bin/
	$(CC) $^ -o $@ -ldl -I./include/ $(ROOT) $(LHAPDF) $(CFLAGS)

clean:
	rm -f bin/*.exe    

