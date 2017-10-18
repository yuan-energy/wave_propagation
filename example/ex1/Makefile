CC = g++ 
CXXFLAGS = -std=c++11 -Wall -g
SRC = \
../../src/WaveField.cpp

INCLUDE_HEADER_LOCATION=\
-I../../src/ \
-I../../utils/

test_deconvolution: test_deconvolution.cpp ${SRC}
	$(CC) -o test_deconvolution.out test_deconvolution.cpp ${CXXFLAGS} ${SRC} ${INCLUDE_HEADER_LOCATION}

clean:
	rm -f test_deconvolution.out
	rm -f wave_at_depth_*.txt