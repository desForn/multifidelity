#pragma once

#include "utility.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Arithmetic
{
    static_assert(std::is_same_v<real_t, float> or std::is_same_v<real_t, double> or
                  std::is_same_v<real_t, long double>, "Arithmetic::real_t: invalid_type");

    extern constexpr real_t default_tolerance =
            std::is_same_v<real_t, float> ? 1E-6F :
            std::is_same_v<real_t, double> ? 1E-14 : 1E-18L;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class type>
    concept real_complex_c = std::is_same_v<std::remove_cvref_t<type>, real_t> or
                             std::is_same_v<std::remove_cvref_t<type>, complex_t>;

    template<class type>
    concept arithmetic_c = std::is_arithmetic_v<type> or real_complex_c<type>;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    constexpr index_t exp2(index_t);
    constexpr integer_t exp2(integer_t);
    constexpr index_t log2(index_t);
    constexpr integer_t log2(integer_t);
    constexpr index_t log2_ceil(index_t);
    constexpr integer_t log2_ceil(integer_t);

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    constexpr auto equal_to(real_t tolerance = default_tolerance);
    constexpr auto not_equal_to(real_t tolerance = default_tolerance);
    constexpr auto greater(real_t tolerance = default_tolerance);
    constexpr auto less(real_t tolerance = default_tolerance);
    constexpr auto greater_equal(real_t tolerance = default_tolerance);
    constexpr auto less_equal(real_t tolerance = default_tolerance);

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    constexpr index_t exp2(index_t arg)
    {
        return 1 << arg;
    }

    constexpr integer_t exp2(integer_t arg)
    {
        ASSERT_ASSUME(arg >= 0);

        return exp2(static_cast<index_t>(arg));
    }

    constexpr index_t log2(index_t arg)
    {
        ASSERT_ASSUME(arg != 0);

        return std::bit_width(arg) - 1;
    }

    constexpr integer_t log2(integer_t arg)
    {
        ASSERT_ASSUME(arg > 0);

        return log2(static_cast<index_t>(arg));
    }

    constexpr index_t log2_ceil(index_t arg)
    {
        ASSERT_ASSUME(arg != 0);

        if (arg == 1)
            return 0;

        return log2(arg - 1) + 1;
    }

    constexpr integer_t log2_ceil(integer_t arg)
    {
        ASSERT_ASSUME(arg > 0);

        if (arg == 1)
            return 0;

        return log2(arg - 1) + 1;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class type0, class type1>
        requires(!real_complex_c<type0> and !real_complex_c<type1>)
        bool equal_to(const type0 &arg0, const type1 &arg1, real_t)
        {
            return arg0 == arg1;
        }

        bool equal_to(real_complex_c auto arg0, real_complex_c auto arg1, real_t tolerance)
        {
            using type0 = decltype(arg0);
            using type1 = decltype(arg1);

            if constexpr (std::is_same_v<type0, real_t> and std::is_same_v<type1, real_t>)
                return std::abs(arg0 - arg1) <= tolerance;

            else
                return std::norm(static_cast<complex_t>(arg0) - static_cast<complex_t>(arg1)) <=
                       tolerance * tolerance;
        }

        template<std::size_t size, real_complex_c type0, real_complex_c type1>
        bool equal_to(const std::array<type0, size> &arg0, const std::array<type1, size> &arg1,
                      real_t tolerance)
        {
            return std::ranges::equal
            (arg0, arg1, [tolerance](type0 i0, type1 i1) -> bool
            {
                return equal_to(i0, i1, tolerance);
            });
        }

        template<index_t i, class... args_types0, class... args_types1>
        bool equal_to_tuple
            (const std::tuple<args_types0...> &arg0, const std::tuple<args_types1...>&arg1,
             real_t tolerance)
        {
            if constexpr (i == sizeof...(args_types0))
                return true;

            else
            {
                using type = std::common_type_t
                        <Utility::get_type<i, args_types0...>,
                         Utility::get_type<i, args_types1...>>;

                if constexpr (real_complex_c<type>)
                {
                    if (!equal_to(std::get<i>(arg0), std::get<i>(arg1), tolerance))
                        return false;
                }

                else
                {
                    if (std::get<i>(arg0) != std::get<i>(arg1))
                        return false;
                }

                return equal_to_tuple<i + 1>(arg0, arg1, tolerance);
            }
        }

        template<class... args_types0, class ...args_types1>
        requires(sizeof...(args_types0) == sizeof...(args_types1))
        bool equal_to(const std::tuple<args_types0...> &arg0, const std::tuple<args_types1...>&arg1,
                      real_t tolerance)
        {
            return equal_to_tuple<0>(arg0, arg1, tolerance);
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

        template<class type0, class type1>
        requires(!real_complex_c<type0> and !real_complex_c<type1>)
        bool less(const type0 &arg0, const type1 &arg1, real_t)
        {
            return arg0 < arg1;
        }

        bool less(real_t arg0, real_t arg1, real_t tolerance)
        {
            return arg0 < arg1 + tolerance;
        }
        
        bool less(complex_t arg0, complex_t arg1, real_t tolerance)
        {
        	if (arg0.real() < arg1.real() + tolerance)
        		return true;
        		
    		else if (arg1.real() < arg0.real() + tolerance)
    			return false;
    			
			else
				return arg0.imag() < arg1.imag() + tolerance;
		}

        template<std::size_t size>
        bool less(const std::array<real_t, size> &arg0, const std::array<real_t, size> &arg1,
                  real_t tolerance)
        {
            for (auto it0 = std::cbegin(arg0), it1 = std::cbegin(arg1); it0 != std::cend(arg0);
                 ++it0, ++it1)
            {
                if (*it0 < *it1 + tolerance)
                    return true;

                if (*it1 < *it0 + tolerance)
                    return false;
            }

            return false;
        }

        template<index_t i, class... args_types0, class... args_types1>
        bool less_tuple
            (const std::tuple<args_types0...> &arg0, const std::tuple<args_types1...>&arg1,
             real_t tolerance)
        {
            if constexpr (i == sizeof...(args_types0))
                return false;

            else
            {
                using type = std::common_type_t
                        <Utility::get_type<i, args_types0...>,
                         Utility::get_type<i, args_types1...>>;

                if constexpr (real_complex_c<type>)
                {
                    if (less(std::get<i>(arg0), std::get<i>(arg1), tolerance))
                        return true;

                    if (!equal_to(std::get<i>(arg0), std::get<i>(arg1), std::abs(tolerance)))
                        return false;
                }

                else
                {
                    if (std::get<i>(arg0) < std::get<i>(arg1))
                        return true;

                    if (std::get<i>(arg0) != std::get<i>(arg1))
                        return false;
                }

                return less_tuple<i + 1>(arg0, arg1, tolerance);
            }
        }

        template<class... args_types0, class... args_types1>
        requires(sizeof...(args_types0) == sizeof...(args_types1))
        bool less(const std::tuple<args_types0...> &arg0, const std::tuple<args_types1...>&arg1,
                  real_t tolerance)
        {
            return less_tuple<0>(arg0, arg1, tolerance);
        }
    }

    constexpr auto equal_to(real_t tolerance)
    {
        return [tolerance](const auto &arg0, const auto &arg1) -> bool
        {
            return Apparatus::equal_to(arg0, arg1, tolerance);
        };
    }

    constexpr auto not_equal_to(real_t tolerance)
    {
        return [tolerance](const auto &arg0, const auto &arg1) -> bool
        {
            return !Apparatus::equal_to(arg0, arg1, tolerance);
        };
    }

    constexpr auto greater(real_t tolerance)
    {
        return [tolerance](const auto &arg0, const auto &arg1) -> bool
        {
            return Apparatus::less(arg1, arg0, tolerance);
        };
    }

    constexpr auto less(real_t tolerance)
    {
        return [tolerance](const auto &arg0, const auto &arg1) -> bool
        {
            return Apparatus::less(arg0, arg1, tolerance);
        };
    }

    constexpr auto greater_equal(real_t tolerance)
    {
        return [tolerance](const auto &arg0, const auto &arg1) -> bool
        {
            return !Apparatus::less(arg0, arg1, tolerance);
        };
    }

    constexpr auto less_equal(real_t tolerance)
    {
        return [tolerance](const auto &arg0, const auto &arg1) -> bool
        {
            return !Apparatus::less(arg1, arg0, tolerance);
        };
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
