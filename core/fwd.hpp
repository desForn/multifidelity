#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <functional>
#include <future>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <ranges>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

#ifdef NDEBUG
#define ASSERT_ASSUME(expression) [[assume(expression)]]
#else
#define ASSERT_ASSUME(expression) assert(expression)

#endif

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
