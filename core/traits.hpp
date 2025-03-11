#pragma once

#include "fwd.hpp"

namespace Traits
{
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class...>
    constexpr bool always_false = false;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class lhs_type, class rhs_type>
    concept comparable_c =
    requires(lhs_type &&lhs, rhs_type &&rhs)
    {
        lhs == rhs;
        rhs == lhs;
        lhs != rhs;
        rhs != lhs;
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    
    namespace Apparatus
    {
        template<class>
        struct span_apparatus
        {
            static constexpr bool value = false;
        };
        
        template<class value_type>
        struct span_apparatus<std::span<value_type>>
        {
            static constexpr bool value = true;
        };
    }
    
    template<class type>
    concept span_c = Apparatus::span_apparatus<type>::value;
    
    template<class container_type>
    concept strict_container_c =
    requires(container_type a, const container_type b)
    {
        requires std::same_as<typename container_type::reference,
                typename container_type::value_type &>;
        requires std::same_as<typename container_type::const_reference,
                const typename container_type::value_type &>;
        requires std::forward_iterator<typename container_type::iterator>;
        requires std::forward_iterator<typename container_type::const_iterator>;
        requires std::signed_integral<typename container_type::difference_type>;
        requires std::unsigned_integral<typename container_type::size_type>;

        requires std::same_as<typename container_type::difference_type,
                typename std::iterator_traits<typename container_type::iterator>::difference_type>;

        requires std::same_as<typename container_type::difference_type,
                typename std::iterator_traits
                        <typename container_type::const_iterator>::difference_type>;

        requires std::regular<container_type>;
        requires std::swappable<container_type>;
        requires std::move_constructible<container_type>;
        requires std::is_move_assignable_v<container_type>;
        requires std::destructible<typename container_type::value_type>;

        { a.begin() } -> std::same_as<typename container_type::iterator>;
        { a.end() } -> std::same_as<typename container_type::iterator>;
        { b.begin() } -> std::same_as<typename container_type::const_iterator>;
        { b.end() } -> std::same_as<typename container_type::const_iterator>;
        { a.cbegin() } -> std::same_as<typename container_type::const_iterator>;
        { a.cend() } -> std::same_as<typename container_type::const_iterator>;
        { b.cbegin() } -> std::same_as<typename container_type::const_iterator>;
        { b.cend() } -> std::same_as<typename container_type::const_iterator>;
        { a.size() } -> std::same_as<typename container_type::size_type>;
        { a.max_size() } -> std::same_as<typename container_type::size_type>;
        { a.empty() } -> std::same_as<bool>;
    };
    
    template<class container_type>
    concept container_c = span_c<container_type> or strict_container_c<container_type>;

    template<class container_type, class value_type>
    concept container_value_type_c = container_c<container_type> and std::is_same_v
            <typename container_type::value_type, value_type>;
            
    template<class range_type, class value_type>
    concept range_c =
        std::ranges::range<range_type> and
        std::same_as<std::ranges::range_value_t<range_type>, value_type>;
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
