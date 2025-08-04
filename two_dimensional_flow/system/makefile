CXX ?= g++
CXXFLAGS ?= -std=c++26 -Wall -Wextra -pedantic -Wno-terminate -O3

all: angle collect_data qoi raw_data smolyak_approximation write_snapshot_data

angle: angle.cpp
	${CXX} angle.cpp -o angle ${CXXFLAGS} -I ../../

collect_data: collect_data.cpp
	${CXX} collect_data.cpp -o collect_data ${CXXFLAGS} -I ../../

qoi: qoi.cpp 
	${CXX} qoi.cpp -o qoi ${CXXFLAGS} -I ../../

raw_data: raw_data.cpp
	${CXX} raw_data.cpp -o raw_data ${CXXFLAGS} -I ../../

smolyak_approximation: smolyak_approximation.cpp
	${CXX} smolyak_approximation.cpp -o smolyak_approximation ${CXXFLAGS} -I ../../

write_snapshot_data: write_snapshot_data.cpp
	${CXX} write_snapshot_data.cpp -o write_snapshot_data ${CXXFLAGS} -I ../../

