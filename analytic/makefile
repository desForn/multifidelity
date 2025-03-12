all: interpolation_1 interpolation_2 interpolation_3 interpolation_4 interpolation_5 quadrature_2 quadrature_4

interpolation_1: interpolation.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 interpolation.cpp -o interpolation_1 -std=c++26 -Wall -Wextra -pedantic -Wno-terminate -O3 -I ../ -D N_VARIATES=1

interpolation_2: interpolation.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 interpolation.cpp -o interpolation_2 -std=c++26 -Wall -Wextra -pedantic -Wno-terminate -O3 -I ../ -D N_VARIATES=2

interpolation_3: interpolation.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 interpolation.cpp -o interpolation_3 -std=c++26 -Wall -Wextra -pedantic -Wno-terminate -O3 -I ../ -D N_VARIATES=3

interpolation_4: interpolation.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 interpolation.cpp -o interpolation_4 -std=c++26 -Wall -Wextra -pedantic -Wno-terminate -O3 -I ../ -D N_VARIATES=4

interpolation_5: interpolation.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 interpolation.cpp -o interpolation_5 -std=c++26 -Wall -Wextra -pedantic -Wno-terminate -O3 -I ../ -D N_VARIATES=5

quadrature_2: quadrature.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 quadrature.cpp -o quadrature_2 -std=c++26 -Wall -Wextra -pedantic -Wno-terminate -O3 -I ../ -D N_VARIATES=1

quadrature_4: quadrature.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 quadrature.cpp -o quadrature_4 -std=c++26 -Wall -Wextra -pedantic -Wno-terminate -O3 -I ../ -D N_VARIATES=3

