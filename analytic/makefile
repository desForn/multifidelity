CXX ?= g++
CXXFLAGS ?= -std=c++26 -Wall -Wextra -pedantic -Wno-terminate -O3

all: interpolation_1 interpolation_2 interpolation_3 interpolation_4 interpolation_5 quadrature_2 quadrature_4

interpolation_1: interpolation.cpp
	${CXX} interpolation.cpp -o interpolation_1 ${CXXFLAGS} -I ../ -D N_VARIATES=1

interpolation_2: interpolation.cpp
	${CXX} interpolation.cpp -o interpolation_2 ${CXXFLAGS} -I ../ -D N_VARIATES=2

interpolation_3: interpolation.cpp
	${CXX} interpolation.cpp -o interpolation_3 ${CXXFLAGS} -I ../ -D N_VARIATES=3

interpolation_4: interpolation.cpp
	${CXX} interpolation.cpp -o interpolation_4 ${CXXFLAGS} -I ../ -D N_VARIATES=4

interpolation_5: interpolation.cpp
	${CXX} interpolation.cpp -o interpolation_5 ${CXXFLAGS} -I ../ -D N_VARIATES=5

quadrature_2: quadrature.cpp
	${CXX} quadrature.cpp -o quadrature_2 ${CXXFLAGS} -I ../ -D N_VARIATES=1

quadrature_4: quadrature.cpp
	${CXX} quadrature.cpp -o quadrature_4 ${CXXFLAGS} -I ../ -D N_VARIATES=3

