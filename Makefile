ROOT=`root-config --cflags` `root-config --libs --glibs` 
#This is the original makefile for when everything was executed from create_templates.cxx

#LHAPDF_DIR=/home/pili/local

#LHAPDF=-I${#LHAPDF_DIR}/include -L${#LHAPDF_DIR}/lib/ -lLHAPDF
#CFLAGS = -O3  -Wall -std=c++14 -pedantic -Wshadow
CFLAGS = -O3  -Wall -std=c++1y -pedantic -Wshadow
#CC=g++
CC=$(CXX)

all: bin/TemplateStruct.o bin/create_templates.exe

bin/create_templates.exe: src/%.cxx, src/%.h
	@mkdir -p tmp/
	@mkdir -p bin/
	$(CC) $^ -o $@ -ldl -I./include/ $(ROOT) $(LHAPDF) $(CFLAGS)

bin/%.o: src/%.cxx src/%.h
	$(CC) $< -o $@ $(ROOT) $(LHAPDF) $(CFLAGS)

clean:
	rm -f bin/*.exe    

