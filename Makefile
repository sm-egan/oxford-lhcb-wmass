ROOT=`root-config --cflags` `root-config --libs --glibs` 

#LHAPDF_DIR=/home/pili/local

#LHAPDF=-I${#LHAPDF_DIR}/include -L${#LHAPDF_DIR}/lib/ -lLHAPDF
#CFLAGS = -O3  -Wall -std=c++14 -pedantic -Wshadow
CFLAGS = -O3  -Wall -std=c++1y -pedantic -Wshadow
#CC=g++
CC=$(CXX)

all: bin/build_compare_histograms.exe  #bin/create_templates.exe

#bin/create_templates.exe: src/create_templates.cxx src/TemplateStruct.cxx src/TemplateStruct.h
#	@mkdir -p tmp/
#	@mkdir -p bin/
#	$(#CC) $^ -o $#@ -ldl -I./include/ $(#ROOT) $(#LHAPDF) $(#CFLAGS)

bin/build_compare_histograms.exe: src/build_compare_histograms.cxx
	@mkdir -p tmp/
	@mkdir -p bin/
	$(CC) $^ -o $@ -ldl -I./include/ $(ROOT) $(LHAPDF) $(CFLAGS)


#bin/%.o: src/%.cxx src/%.h
#	$(#CC) $< -o $#@ $(#ROOT) $(#LHAPDF) $(#CFLAGS)
#bin/Zsample_analysis.exe: src/Zsample_analysis.cxx src/build_compare_histograms.cxx src/build_compare_histograms.h
#	@mkdir -p tmp/
#	@mkdir -p bin/
#	$(#CC) $^ -o $@ -ldl -I./include/ $(#ROOT) $(#LHAPDF) $(#CFLAGS)

clean:
	rm -f bin/*.exe    

