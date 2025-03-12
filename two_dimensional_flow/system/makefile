all: angle collect_data qoi raw_data smolyak_approximation write_snapshot_data

angle: angle.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 angle.cpp -o angle -std=c++26 -Wall -Wextra -pedantic -Wno-terminate -O3 -I ../../

collect_data: collect_data.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 collect_data.cpp -o collect_data -std=c++26 -Wall -Wextra -pedantic -Wno-terminate -O3 -I ../../

qoi: qoi.cpp 
	/usr/local/gcc/gcc-14/bin/g++-14 qoi.cpp -o qoi -std=c++26 -Wall -Wextra -pedantic -Wno-terminate -O3 -I ../../

raw_data: raw_data.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 raw_data.cpp -o raw_data -std=c++26 -Wall -Wextra -pedantic -Wno-terminate -O3 -I ../../

smolyak_approximation: smolyak_approximation.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 smolyak_approximation.cpp -o smolyak_approximation -std=c++26 -Wall -Wextra -pedantic -Wno-terminate -O3 -I ../../

write_snapshot_data: write_snapshot_data.cpp
	/usr/local/gcc/gcc-14/bin/g++-14 write_snapshot_data.cpp -o write_snapshot_data -std=c++26 -Wall -Wextra -pedantic -Wno-terminate -O3 -I ../../

