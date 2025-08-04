CXX ?= g++
CXXFLAGS ?= -std=c++26 -Wall -Wextra -pedantic -O3

all: collect_data generate_bcs generate_nu qoi platoon_truncation platoon_richardson_extrapolation time

collect_data: collect_data.cpp
	${CXX} collect_data.cpp -o collect_data ${CXXFLAGS} -I ../../

generate_bcs: generate_bcs.cpp
	${CXX} generate_bcs.cpp -o generate_bcs ${CXXFLAGS} -I ../../

generate_nu: generate_nu.cpp
	${CXX} generate_nu.cpp -o generate_nu ${CXXFLAGS} -I ../../

qoi: qoi.cpp
	${CXX} qoi.cpp -o qoi ${CXXFLAGS} -I ../../

platoon_truncation: platoon_truncation.cpp
	${CXX} platoon_truncation.cpp -o platoon_truncation ${CXXFLAGS} -I ../../

platoon_richardson_extrapolation: platoon_richardson_extrapolation.cpp
	${CXX} platoon_richardson_extrapolation.cpp -o platoon_richardson_extrapolation ${CXXFLAGS} -I ../../

time:
	${CXX} time.cpp -o time -Ofast -fno-finite-math-only -I ../../ -DNDEBUG -std=c++26 -march=native -flto

