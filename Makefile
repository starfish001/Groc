all:  gretc
gretc: gzstream.cpp gretc.hash.mul_thread.cpp
	g++   -O3 -static  -o $@ $^ -lz -pthread
