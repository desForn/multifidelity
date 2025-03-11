#pragma once

#include "arithmetic.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Utility
{
    template<Arithmetic::index_t n_variates>
    class counter
    {
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using signed_indices_type = std::array<Arithmetic::integer_t, n_variates>;
        using indices_type = std::array<Arithmetic::index_t, n_variates>;

    private:
        struct counter_iterator;
        struct void_struct {};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        indices_type begin_{};
        indices_type end_;
        signed_indices_type strides_{uniform_array<Arithmetic::integer_t, n_variates>(1)};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        counter() = delete;
        ~counter() = default;

        counter(const counter &) = delete;
        counter &operator=(const counter &) = delete;

        counter(counter &&) noexcept = default;
        counter &operator=(counter &&) noexcept = default;

        counter(Arithmetic::index_t end);
        counter(const indices_type &end);
        counter(const indices_type &begin, const indices_type &end);
        counter(const indices_type &begin, const indices_type &end, const indices_type &strides);
        counter(const indices_type &begin, const indices_type &end,
                const signed_indices_type &strides);

    public:
        counter_iterator begin() const;
        void_struct end() const;
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Arithmetic::index_t n_variates>
    struct counter<n_variates>::counter_iterator
    {
        indices_type count_{};
        const indices_type &begin_;
        const indices_type &end_;
        const signed_indices_type &strides_;

        counter_iterator
                (const indices_type &begin, const indices_type &end,
                 const signed_indices_type &strides);

        bool operator!=(void_struct);
        const indices_type &operator*();
        counter_iterator &operator++();
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Arithmetic::index_t n_variates>
    counter<n_variates>::counter(Arithmetic::index_t end) :
        counter{Utility::uniform_array<Arithmetic::index_t, n_variates>(end)} {}

    template<Arithmetic::index_t n_variates>
    counter<n_variates>::counter(const indices_type &end) : end_{end} {}

    template<Arithmetic::index_t n_variates>
    counter<n_variates>::counter(const indices_type &begin, const indices_type &end) :
            begin_{begin}, end_{end}
    {
        ASSERT_ASSUME(std::ranges::equal(begin_, end_, [](auto i, auto j) { return i <= j; }));

        return;
    }

    template<Arithmetic::index_t n_variates>
    counter<n_variates>::counter
            (const indices_type &begin, const indices_type &end, const indices_type &strides) :
            counter(begin, end, Utility::signed_array(strides)) {}

    template<Arithmetic::index_t n_variates>
    counter<n_variates>::counter
            (const indices_type &begin, const indices_type &end,
             const signed_indices_type &strides) :
            begin_{begin}, end_{end}, strides_{strides}
    {
        ASSERT_ASSUME(std::ranges::find(strides_, 0) == std::cend(strides_));

        for (Arithmetic::index_t i = 0; i != n_variates; ++i)
            ASSERT_ASSUME(strides_[i] *
                          static_cast<Arithmetic::integer_t>(end_[i] - begin_[i]) >= 0);

        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Arithmetic::index_t n_variates>
    auto counter<n_variates>::begin() const -> counter_iterator
    {
        return counter_iterator{begin_, end_, strides_};
    }

    template<Arithmetic::index_t n_variates>
    auto counter<n_variates>::end() const -> void_struct
    {
        return void_struct{};
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<Arithmetic::index_t n_variates>
    counter<n_variates>::counter_iterator::counter_iterator
            (const indices_type &begin, const indices_type &end,
             const signed_indices_type &strides) :
            count_{begin}, begin_{begin}, end_{end}, strides_{strides} {}

    template<Arithmetic::index_t n_variates>
    bool counter<n_variates>::counter_iterator::operator!=(void_struct)
    {
        for (Arithmetic::index_t i = 0; i != n_variates; ++i)
            if (strides_[i] * static_cast<Arithmetic::integer_t>(end_[i] - count_[i]) <= 0)
                return false;

        return true;
    }

    template<Arithmetic::index_t n_variates>
    auto counter<n_variates>::counter_iterator::operator*() -> const indices_type &
    {
        return count_;
    }

    template<Arithmetic::index_t n_variates>
    auto counter<n_variates>::counter_iterator::operator++() -> counter_iterator &
    {
        count_.back() += strides_.back();

        for (Arithmetic::index_t i = n_variates - 1; i != 0; --i)
        {
            if (strides_[i] * static_cast<Arithmetic::integer_t>(end_[i] - count_[i]) <= 0)
            {
                count_[i] = begin_[i];
                count_[i - 1] += strides_[i - 1];
            }

            else
                break;
        }

        return *this;
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
