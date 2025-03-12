#pragma once

#include "chebyshev_polynomial.hpp"
#include "multivariate_polynomial_point_value.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Polynomial
{
    void discrete_fourier_transform(std::span<complex_t>);
    void inverse_discrete_fourier_transform(std::span<complex_t>);

    void multidimensional_discrete_fourier_transform
        (std::span<complex_t>, std::span<const index_t>);
    void inverse_multidimensional_discrete_fourier_transform
        (std::span<complex_t>, std::span<const index_t>);

    chebyshev_polynomial_i_kind discrete_cosine_transform(const chebyshev_polynomial_point_value &);
    chebyshev_polynomial_point_value discrete_cosine_transform(const chebyshev_polynomial_i_kind &);

    template<index_t n_variates>
    multivariate_polynomial<chebyshev_polynomial_i_kind, n_variates> discrete_cosine_transform
            (const multivariate_polynomial<chebyshev_polynomial_point_value, n_variates> &);

    template<index_t n_variates>
    multivariate_polynomial<chebyshev_polynomial_point_value, n_variates> discrete_cosine_transform
            (const multivariate_polynomial<chebyshev_polynomial_i_kind, n_variates> &);

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        std::uint8_t reverse_bit(std::uint8_t n)
        {
            n = (n & 0xF0) >> 4 | (n & 0x0F) << 4;
            n = (n & 0xCC) >> 2 | (n & 0x33) << 2;
            n = (n & 0xAA) >> 1 | (n & 0x55) << 1;

            return n;
        }

        std::uint16_t reverse_bit(std::uint16_t n)
        {
            n = (n & 0xFF00) >> 8 | (n & 0x00FF) << 8;
            n = (n & 0xF0F0) >> 4 | (n & 0X0F0F) << 4;
            n = (n & 0xCCCC) >> 2 | (n & 0x3333) << 2;
            n = (n & 0xAAAA) >> 1 | (n & 0X5555) << 1;

            return n;
        }

        std::uint32_t reverse_bit(std::uint32_t n)
        {
            n = (n & 0xFFFF0000) >> 16 | (n & 0x0000FFFF) << 16;
            n = (n & 0xFF00FF00) >> 8  | (n & 0x00FF00FF) << 8;
            n = (n & 0xF0F0F0F0) >> 4  | (n & 0X0F0F0F0F) << 4;
            n = (n & 0xCCCCCCCC) >> 2  | (n & 0x33333333) << 2;
            n = (n & 0xAAAAAAAA) >> 1  | (n & 0X55555555) << 1;

            return n;
        }

        std::uint64_t reverse_bit(std::uint64_t n)
        {
            n = (n & 0xFFFFFFFF00000000) >> 32 | (n & 0x00000000FFFFFFFF) << 32;
            n = (n & 0xFFFF0000FFFF0000) >> 16 | (n & 0x0000FFFF0000FFFF) << 16;
            n = (n & 0xFF00FF00FF00FF00) >> 8  | (n & 0x00FF00FF00FF00FF) << 8;
            n = (n & 0xF0F0F0F0F0F0F0F0) >> 4  | (n & 0X0F0F0F0F0F0F0F0F) << 4;
            n = (n & 0xCCCCCCCCCCCCCCCC) >> 2  | (n & 0x3333333333333333) << 2;
            n = (n & 0xAAAAAAAAAAAAAAAA) >> 1  | (n & 0X5555555555555555) << 1;

            return n;
        }

        index_t reverse_bit(index_t n, index_t b)
        {
            return reverse_bit(static_cast<uint64_t>(n)) >> (64 - b);
        }
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    void discrete_fourier_transform(std::span<complex_t> ret)
    {
        index_t n = std::size(ret);
        ASSERT_ASSUME(n != 0 and (n & (n - 1)) == 0); // Assert n is arg power of 2

        index_t log_n = Arithmetic::log2(n);

        for (index_t i = 0; i != n; ++i)
            if (index_t j = Apparatus::reverse_bit(i, log_n); j > i)
                std::swap(ret[i], ret[j]);

        for (index_t i = 1; i != log_n + 1; ++i)
        {
            index_t m = Arithmetic::exp2(i);

            real_t alpha = static_cast<real_t>(-2) * static_cast<real_t>(M_PI) /
                           static_cast<real_t>(m);

            complex_t wm{std::cos(alpha), std::sin(alpha)};

            for (index_t j = 0; j != n; j += m)
            {
                complex_t w{1};

                for (index_t k = 0; k != m / 2; ++k)
                {
                    complex_t t = w * ret[j + k + m / 2];
                    complex_t u = ret[j + k];
                    ret[j + k] = u + t;
                    ret[j + k + m / 2] = u - t;
                    w *= wm;
                }
            }
        }

        std::ranges::for_each(ret, [n](complex_t &i) { i /= n; });
    }

    void inverse_discrete_fourier_transform(std::span<complex_t> ret)
    {
        index_t n = std::size(ret);
        ASSERT_ASSUME(n != 0 and (n & (n - 1)) == 0); // Assert n is arg power of 2

        index_t log_n = Arithmetic::log2(n);

        for (index_t i = 1; i < n / 2; ++i)
            std::swap(ret[i], ret[n - i]);

        for (index_t i = 0; i != n; ++i)
            if (index_t j = Apparatus::reverse_bit(i, log_n); j > i)
                std::swap(ret[i], ret[j]);

        for (index_t i = 1; i != log_n + 1; ++i)
        {
            index_t m = 1 << i;

            real_t alpha = static_cast<real_t>(-2) * static_cast<real_t>(M_PI) /
                           static_cast<real_t>(m);

            complex_t wm{std::cos(alpha), std::sin(alpha)};

            for (index_t j = 0; j != n; j += m)
            {
                complex_t w{1};

                for (index_t k = 0; k != m / 2; ++k)
                {
                    complex_t t = w * ret[j + k + m / 2];
                    complex_t u = ret[j + k];
                    ret[j + k] = u + t;
                    ret[j + k + m / 2] = u - t;
                    w *= wm;
                }
            }
        }
    }

    void multidimensional_discrete_fourier_transform
            (std::span<complex_t> arg, std::span<const index_t> dimensions)
    {
        ASSERT_ASSUME(!arg.empty() and !dimensions.empty() and
                      std::size(arg) == std::accumulate(
                              std::cbegin(dimensions), std::cend(dimensions),
                              static_cast<index_t>(1), std::multiplies<index_t>{}));

        if (std::size(dimensions) == 1)
        {
            discrete_fourier_transform(arg);
            return;
        }

        index_t stride = std::size(arg) / dimensions.front();

        std::span<const index_t> tail_dimensions
                {std::cbegin(dimensions) + 1, std::cend(dimensions)};

        // DFT of each row
        for (auto it = std::begin(arg); it != std::end(arg); it += stride)
            multidimensional_discrete_fourier_transform
                (std::span<complex_t>{it, it + stride}, tail_dimensions);

        // v = transpose(arg)
        std::vector<complex_t> v(std::size(arg));

        auto it_0 = std::begin(arg);
        auto it_1 = std::begin(v);

        for (; it_0 != std::begin(arg) + stride; ++it_0, it_1 += dimensions.front())
        {
            auto it_0_bis = it_0;
            auto it_1_bis = it_1;

            /* The loop should read as follows:
             *  for (; it_0_bis < std::end(arg); it_0_bis += stride, ++it_1_bis)
             *      *it_1_bis = *it_0_bis;
             *
             * However, it_0_bis is incremented past std::end(arg).
             * Although it is never dereferenced, wit GLIBCXX_DEBUG mode it gives an error. */

            for (; it_0_bis < std::end(arg) - stride; it_0_bis += stride, ++it_1_bis)
                *it_1_bis = *it_0_bis;

            *it_1_bis = *it_0_bis;

            // DFT of each row of v, i.e. each column of arg
            discrete_fourier_transform(std::span<complex_t>{it_1, it_1 + dimensions.front()});
        }

        // arg = transpose(v)
        it_0 = std::begin(arg);
        it_1 = std::begin(v);

        for (; it_1 != std::begin(v) + dimensions.front(); ++it_1, it_0 += stride)
        {
            auto it_0_bis = it_0;
            auto it_1_bis = it_1;

            /* The loop should read as follows:
             *  for (; it_1_bis < std::end(v); it_1_bis += dimensions.front(), ++it_0_bis)
             *      *it_0_bis = *it_1_bis;
             *
             * However, it_1_bis is incremented past std::end(v).
             * Although it is never dereferenced, wit GLIBCXX_DEBUG mode it gives an error. */

            for (; it_1_bis < std::end(v) - dimensions.front();
                   it_1_bis += dimensions.front(), ++it_0_bis)
                *it_0_bis = *it_1_bis;

            *it_0_bis = *it_1_bis;
        }

        return;
    }

    void inverse_multidimensional_discrete_fourier_transform
            (std::span<complex_t> arg, std::span<const index_t> dimensions)
    {
        ASSERT_ASSUME(!arg.empty() and !dimensions.empty() and
                      std::size(arg) == std::accumulate(
                              std::cbegin(dimensions), std::cend(dimensions),
                              static_cast<index_t>(1), std::multiplies<index_t>{}));

        if (std::size(dimensions) == 1)
        {
            inverse_discrete_fourier_transform(arg);
            return;
        }

        index_t stride = std::size(arg) / dimensions.front();

        std::span<const index_t> tail_dimensions
            {std::cbegin(dimensions) + 1, std::cend(dimensions)};

        // DFT of each row
        for (auto it = std::begin(arg); it != std::end(arg); it += stride)
            inverse_multidimensional_discrete_fourier_transform
                (std::span<complex_t>{it, it + stride}, tail_dimensions);

        // v = transpose(arg)
        std::vector<complex_t> v(std::size(arg));

        auto it_0 = std::begin(arg);
        auto it_1 = std::begin(v);

        for (; it_0 != std::begin(arg) + stride; ++it_0, it_1 += dimensions.front())
        {
            auto it_0_bis = it_0;
            auto it_1_bis = it_1;

            /* The loop should read as follows:
             *  for (; it_0_bis < std::end(arg); it_0_bis += stride, ++it_1_bis)
             *      *it_1_bis = *it_0_bis;
             *
             * However, it_0_bis is incremented past std::end(arg).
             * Although it is never dereferenced, wit GLIBCXX_DEBUG mode it gives an error. */

            for (; it_0_bis < std::end(arg) - stride; it_0_bis += stride, ++it_1_bis)
                *it_1_bis = *it_0_bis;

            *it_1_bis = *it_0_bis;

            // DFT of each row of v, i.e. each column of arg
            inverse_discrete_fourier_transform
                (std::span<complex_t>{it_1, it_1 + dimensions.front()});
        }

        // arg = transpose(v)
        it_0 = std::begin(arg);
        it_1 = std::begin(v);

        for (; it_1 != std::begin(v) + dimensions.front(); ++it_1, it_0 += stride)
        {
            auto it_0_bis = it_0;
            auto it_1_bis = it_1;

            /* The loop should read as follows:
             *  for (; it_1_bis < std::end(v); it_1_bis += dimensions.front(), ++it_0_bis)
             *      *it_0_bis = *it_1_bis;
             *
             * However, it_1_bis is incremented past std::end(v).
             * Although it is never dereferenced, wit GLIBCXX_DEBUG mode it gives an error. */

            for (; it_1_bis < std::end(v) - dimensions.front();
                   it_1_bis += dimensions.front(), ++it_0_bis)
                *it_0_bis = *it_1_bis;

            *it_0_bis = *it_1_bis;
        }

        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    chebyshev_polynomial_i_kind discrete_cosine_transform
        (const chebyshev_polynomial_point_value &arg)
    {
        index_t extents_arg{arg.extents()};

        if (extents_arg == 0)
            return chebyshev_polynomial_i_kind{};

        if (extents_arg == 1)
            return chebyshev_polynomial_i_kind{arg.f().front()};

        std::vector<complex_t> unfolded(2 * (extents_arg - 1));

        for (index_t i = 0; i != std::size(unfolded); ++i)
            unfolded[i] = arg.f(std::abs(static_cast<integer_t>(extents_arg) - 1 -
                                         static_cast<integer_t>(i)));

        discrete_fourier_transform(std::span{unfolded});

        chebyshev_polynomial_i_kind ret;
        ret.set_order(arg.order());

        auto it_f = std::cbegin(unfolded);
        auto it_ret = std::cbegin(ret.coefficients());

        *it_ret = it_f->real();

        for (++it_ret, ++it_f; it_ret != std::cend(ret.coefficients()) - 1; ++it_ret, ++it_f)
            *it_ret = 2 * it_f->real();

        *it_ret = it_f->real() * (1 + (ret.extents() != extents_arg));

        return ret;
    }

    chebyshev_polynomial_point_value discrete_cosine_transform
            (const chebyshev_polynomial_i_kind &arg)
    {
        index_t extents_arg = arg.extents();

        if (extents_arg == 0)
            return chebyshev_polynomial_point_value{};

        if (extents_arg == 1)
            return chebyshev_polynomial_point_value{arg.coefficients().front()};

        index_t sampling_level =
                Sampling::required_level_exponential_sampling(extents_arg);
        index_t extents_ret = Sampling::n_points_exponential_sampling(sampling_level);

        std::vector<complex_t> vector_unfolded(2 * (extents_ret - 1));
        vector_unfolded[0] = arg[0];

        for (index_t i = 1; i != extents_arg; ++i)
        {
            vector_unfolded[i] = arg[i] / 2;
            vector_unfolded[2 * (extents_arg - 1) - i] = arg[i] / 2;
        }

        if (extents_arg == extents_ret)
            vector_unfolded[extents_ret - 1] = arg[extents_ret - 1];

        inverse_discrete_fourier_transform(std::span{vector_unfolded});

        std::vector<real_t> v(extents_ret);

        for (index_t i = 0; i != extents_ret; ++i)
            v[extents_ret - 1 - i] = vector_unfolded[i].real();

        return chebyshev_polynomial_point_value{std::move(v), sampling_level, arg.order()};
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<index_t n_variates>
    multivariate_polynomial<chebyshev_polynomial_i_kind, n_variates> discrete_cosine_transform
            (const multivariate_polynomial<chebyshev_polynomial_point_value, n_variates> &arg)
    {
        using ret_type = multivariate_polynomial<chebyshev_polynomial_i_kind, n_variates>;
        using indices_type = std::array<index_t, n_variates>;

        index_t size_arg = std::size(arg.f());
        const indices_type &extents_arg = arg.extents();

        if (size_arg == 0)
            return ret_type{};

        if (size_arg == 1)
            return ret_type{arg.f().front()};

        indices_type unfolded_extents = Utility::transform_array(
                [](index_t i) -> index_t
                {
                    return 2 * (i - 1) + (i == 1);
                },
                extents_arg);

        indices_type unfolded_strides;
        unfolded_strides.back() = 1;

        for (index_t i = n_variates - 1; i != 0; --i)
            unfolded_strides[i - 1] = unfolded_strides[i] * unfolded_extents[i];

        std::vector<complex_t> unfolded(unfolded_extents.front() * unfolded_strides.front());

        auto it_unfolded = std::begin(unfolded);

        for (const auto &unfolded_index: Utility::counter(unfolded_extents))
        {
            indices_type folded_index;

            std::transform(std::cbegin(unfolded_index), std::cend(unfolded_index),
                           std::cbegin(arg.extents()),
                           std::begin(folded_index),
                           [](integer_t j, integer_t extent) -> index_t
                           {
                               return std::abs((extent - 1) - j);
                           });

            *it_unfolded++ = arg.f(folded_index);
        }

        multidimensional_discrete_fourier_transform(std::span<complex_t>{unfolded},
                                                    std::span<const index_t>{unfolded_extents});

        ret_type ret;
        ret.set_order(arg.order());

        auto it_ret = std::begin(ret.coefficients());

        for (const auto &ret_index: Utility::counter(ret.extents()))
        {
            index_t unfolded_offset = std::inner_product
                    (std::cbegin(unfolded_strides), std::cend(unfolded_strides),
                     std::cbegin(ret_index), static_cast<index_t>(0));

            real_t c = unfolded[unfolded_offset].real();

            for (index_t v = 0; v != n_variates; ++v)
                if (ret_index[v] != 0 and ret_index[v] != ret.extents(v) - 1)
                    c *= 2;

            *it_ret++ = c;
        }

        return ret;
    }

    template<index_t n_variates>
    multivariate_polynomial<chebyshev_polynomial_point_value, n_variates> discrete_cosine_transform
            (const multivariate_polynomial<chebyshev_polynomial_i_kind, n_variates> &arg)
    {
        using ret_type = multivariate_polynomial<chebyshev_polynomial_point_value, n_variates>;
        using indices_type = std::array<index_t, n_variates>;

        if (std::size(arg.coefficients()) == 0)
            return ret_type{};

        if (std::size(arg.coefficients()) == 1)
            return ret_type{arg.coefficients().front()};

        const indices_type &extents_arg = arg.extents();

        indices_type sampling_level = Utility::transform_array(
                [](index_t i) -> index_t
                {
                    return Sampling::required_level_exponential_sampling(i);
                },
                extents_arg);

        indices_type extents_ret = Utility::transform_array(
                [](index_t i) -> index_t
                {
                    return Sampling::n_points_exponential_sampling(i);
                },
                sampling_level);

        indices_type strides_ret;
        strides_ret.back() = 1;

        for (index_t i = n_variates - 1; i != 0; --i)
            strides_ret[i - 1] = strides_ret[i] * extents_ret[i];

        indices_type extents_unfolded = Utility::transform_array(
                [](index_t i) -> index_t
                {
                    return 2 * (i - 1) + (i == 1);
                },
                extents_ret);

        indices_type strides_unfolded;
        strides_unfolded.back() = 1;

        for (index_t i = n_variates - 1; i != 0; --i)
            strides_unfolded[i - 1] = strides_unfolded[i] * extents_unfolded[i];

        std::vector<complex_t> vector_unfolded(extents_unfolded.front() * strides_unfolded.front());

        auto it_unfolded = std::begin(vector_unfolded);

        for (const auto &unfolded_index: Utility::counter(extents_unfolded))
        {
            indices_type folded_index;

            std::transform(std::cbegin(unfolded_index), std::cend(unfolded_index),
                           std::cbegin(extents_ret),
                           std::begin(folded_index),
                           [](integer_t j, integer_t extent) -> index_t
                           {
                               return extent - 1 - std::abs((extent - 1) - j);
                           });

            if (!std::ranges::equal(folded_index, extents_arg,
                                    [](index_t i, index_t j) -> bool
                                    {
                                        return i < j;
                                    }))
                continue;

            real_t c = arg[folded_index];

            for (index_t v = 0; v != n_variates; ++v)
                if (folded_index[v] != 0 and folded_index[v] != extents_ret[v] - 1)
                    c /= 2;

            *it_unfolded++ = c;
        }

        inverse_multidimensional_discrete_fourier_transform
                (std::span<complex_t>{vector_unfolded}, std::span<const index_t>{extents_unfolded});

        std::vector<real_t> vector_ret(extents_ret.front() * strides_ret.front());

        auto it_ret = std::begin(vector_ret);

        for (const auto &index_ret: Utility::counter(extents_ret))
        {
            indices_type index_mirrored;

            std::transform(std::cbegin(extents_ret), std::cend(extents_ret),
                           std::cbegin(index_ret), std::begin(index_mirrored),
                           [](index_t e, index_t i) -> index_t
                           {
                               return e - 1 - i;
                           });

            index_t offset_unfolded = std::inner_product
                    (std::cbegin(strides_unfolded), std::cend(strides_unfolded),
                     std::cbegin(index_mirrored), static_cast<index_t>(0));

            *it_ret++ = vector_unfolded[offset_unfolded].real();
        }

        ret_type ret{std::move(vector_ret), sampling_level, arg.order()};

        return ret;
    }
}
