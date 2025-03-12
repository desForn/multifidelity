all: collect_data generate_bcs generate_nu qoi platoon_truncation platoon_richardson_extrapolation

collect_data: collect_data.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 collect_data.cpp -o collect_data -std=c++26 -Wall -Wextra -pedantic -O3 -I ../../

generate_bcs: generate_bcs.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 generate_bcs.cpp -o generate_bcs -std=c++26 -Wall -Wextra -pedantic -O3 -I ../../

generate_nu: generate_nu.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 generate_nu.cpp -o generate_nu -std=c++26 -Wall -Wextra -pedantic -O3 -I ../../

qoi: qoi.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 qoi.cpp -o qoi -std=c++26 -Wall -Wextra -pedantic -O3 -I ../../

platoon_truncation: platoon_truncation.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 platoon_truncation.cpp -o platoon_truncation -std=c++26 -Wall -Wextra -pedantic -O3 -I ../../

platoon_richardson_extrapolation: platoon_richardson_extrapolation.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 platoon_richardson_extrapolation.cpp -o platoon_richardson_extrapolation -std=c++26 -Wall -Wextra -pedantic -O3 -I ../../

