#pragma once

#include "fwd.hpp"
#include "invocable/invocable_function_hash_table.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Sampling
{
    template<sampling_traits_c sampling_traits_, index_t offset_ = 0>
    class sampling
    {
    public:
        using sampling_traits = sampling_traits_;

        using coordinates_type = sampling_traits::coordinates_type;
        static constexpr bool new_level_information = sampling_traits_new_level_c<sampling_traits>;
        static constexpr index_t offset = offset_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        class factory_apparatus;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        const std::vector<coordinates_type> coordinates_;
        const std::vector<coordinates_type> new_coordinates_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        sampling() = delete;
        ~sampling() = default;

        sampling(const sampling &) = default;
        sampling &operator=(const sampling &) = default;

        sampling(sampling &&) noexcept = default;
        sampling &operator=(sampling &&) noexcept = default;

        explicit sampling(index_t);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        static const sampling &factory(index_t);
        static void clear_factory();

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        static constexpr index_t n_points(index_t);
        static constexpr index_t n_new_points(index_t);
        static constexpr index_t required_level(index_t);

        static index_t index(const coordinates_type &);
        static coordinates_type coordinate_from_index(index_t);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        [[nodiscard]] index_t level() const;
        [[nodiscard]] index_t n_points() const;
        [[nodiscard]] index_t n_new_points() const;
        [[nodiscard]] std::span<const coordinates_type> coordinates() const &;
        [[nodiscard]] std::vector<coordinates_type> &&coordinates() &&;
        [[nodiscard]] const coordinates_type &operator[](index_t) const;
        [[nodiscard]] std::span<const coordinates_type> new_coordinates() const &;
        [[nodiscard]] std::vector<coordinates_type> &&new_coordinates() &&;

    private:
        [[nodiscard]] std::vector<coordinates_type> new_coordinates_initialiser() const;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        auto begin() const { return std::begin(coordinates_); }
        auto end() const { return std::end(coordinates_); }

        auto cbegin() const { return std::cbegin(coordinates_); }
        auto cend() const { return std::cend(coordinates_); }

        auto rbegin() const { return std::rbegin(coordinates_); }
        auto rend() const { return std::rend(coordinates_); }

        auto crbegin() const { return std::crbegin(coordinates_); }
        auto crend() const { return std::crend(coordinates_); }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<sampling_traits_c sampling_traits, index_t offset>
    class sampling<sampling_traits, offset>::factory_apparatus
    {
        friend class sampling;

        inline static auto lambda_ = [](index_t level)
        {
            return std::make_unique<sampling>(sampling{level});
        };

        inline static Invocable::invocable_function_hash_table invocable_
                {lambda_,
                 Invocable::invocable<const std::unique_ptr<sampling> &(index_t)>{}};
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<sampling_traits_c sampling_traits, index_t offset>
    sampling<sampling_traits, offset>::sampling(index_t level) :
    coordinates_{sampling_traits::coordinates(level + offset)},
    new_coordinates_{this->new_coordinates_initialiser()} {}

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<sampling_traits_c sampling_traits, index_t offset>
    auto sampling<sampling_traits, offset>::factory(index_t level) -> const sampling &
    {
        return *factory_apparatus::invocable_(level);
    }

    template<sampling_traits_c sampling_traits, index_t offset>
    void sampling<sampling_traits, offset>::clear_factory()
    {
        factory_apparatus::invocable_.clear();
        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<sampling_traits_c sampling_traits, index_t offset>
    constexpr index_t sampling<sampling_traits, offset>::n_points(index_t level)
    {
        return sampling_traits::n_points(level + offset);
    }

    template<sampling_traits_c sampling_traits, index_t offset>
    constexpr index_t sampling<sampling_traits, offset>::n_new_points(index_t level)
    {
        level += offset;

        if constexpr (new_level_information)
            return sampling_traits::n_new_points(level);

        else
        {
            auto less = Arithmetic::less(-Arithmetic::default_tolerance);
            std::set<coordinates_type, decltype(less)> set_{less};

            for (index_t level_secundum = 0; level_secundum != level; ++level_secundum)
            {
                const sampling &sampling_ = sampling::factory(level_secundum);

                for (const coordinates_type &i : sampling_)
                    set_.insert(i);
            }

            index_t old_size = std::size(set_);

            const sampling &sampling_ = sampling::factory(level);

            for (const coordinates_type &i : sampling_)
                set_.insert(i);

            return std::size(set_) - old_size;
        }
    }

    template<sampling_traits_c sampling_traits, index_t offset>
    constexpr index_t sampling<sampling_traits, offset>::required_level(index_t n_points)
    {
        index_t ret = sampling_traits::required_level(n_points);

        if (ret < offset)
            return 0;

        return ret - offset;
    }

    template<sampling_traits_c sampling_traits, index_t offset>
    index_t sampling<sampling_traits, offset>::index(const coordinates_type &arg)
    {
        auto equal_to = Arithmetic::equal_to();

        if constexpr (new_level_information)
        {
            index_t ret = 0;

            for (index_t level = 0; true; ++level)
            {
                const sampling &sampling_ = sampling::factory(level);

                for (index_t i = 0; i != sampling_.n_points(); ++i)
                    if (sampling_traits::is_new_point(i, level))
                    {
                        if (equal_to(arg, sampling_[i]))
                            return ret;

                        ++ret;
                    }
            }
        }

        else
        {
            auto less = Arithmetic::less(-Arithmetic::default_tolerance);
            std::set<coordinates_type, decltype(less)> set_{less};

            for (index_t level = 0; true; ++level)
            {
                const sampling &sampling_ = sampling::factory(level);

                for (const coordinates_type &i : sampling_)
                {
                    if (equal_to(arg, i))
                        return std::size(set_);

                    set_.insert(i);
                }
            }
        }
    }

    template<sampling_traits_c sampling_traits, index_t offset>
    auto sampling<sampling_traits, offset>::coordinate_from_index(index_t arg) -> coordinates_type
    {
        if constexpr (new_level_information)
        {
            index_t level = 0;
            for (; n_new_points(level) <= arg; arg -= n_new_points(level++)) ;

            return sampling::factory(level).new_coordinates()[arg];
        }

        else
        {
            ++arg;

            auto less = Arithmetic::less(-Arithmetic::default_tolerance);
            std::set<coordinates_type, decltype(less)> set_{less};

            for (index_t level = 0; true; ++level)
            {
                const sampling &sampling_ = sampling::factory(level);

                for (const coordinates_type &i : sampling_)
                {
                    set_.insert(i);

                    if (std::size(set_) == arg)
                        return i;
                }
            }
        }
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<sampling_traits_c sampling_traits, index_t offset>
    index_t sampling<sampling_traits, offset>::level() const
    {
        return required_level(n_points());
    }

    template<sampling_traits_c sampling_traits, index_t offset>
    index_t sampling<sampling_traits, offset>::n_points() const
    {
        return std::size(coordinates_);
    }

    template<sampling_traits_c sampling_traits, index_t offset>
    index_t sampling<sampling_traits, offset>::n_new_points() const
    {
        return n_new_points(level());
    }

    template<sampling_traits_c sampling_traits, index_t offset>
    auto sampling<sampling_traits, offset>::coordinates() const &
    -> std::span<const coordinates_type>
    {
        return {coordinates_};
    }

    template<sampling_traits_c sampling_traits, index_t offset>
    auto sampling<sampling_traits, offset>::coordinates() && -> std::vector<coordinates_type> &&
    {
        return std::move(coordinates_);
    }

    template<sampling_traits_c sampling_traits, index_t offset>
    auto sampling<sampling_traits, offset>::operator[](index_t arg) const
    -> const coordinates_type &
    {
        return coordinates_[arg];
    }

    template<sampling_traits_c sampling_traits, index_t offset>
    auto sampling<sampling_traits, offset>::new_coordinates() const &
    -> std::span<const coordinates_type>
    {
        return {new_coordinates_};
    }

    template<sampling_traits_c sampling_traits, index_t offset>
    auto sampling<sampling_traits, offset>::new_coordinates() && -> std::vector<coordinates_type> &&
    {
        return std::move(new_coordinates_);
    }

    template<sampling_traits_c sampling_traits, index_t offset>
    auto sampling<sampling_traits, offset>::new_coordinates_initialiser() const
    -> std::vector<coordinates_type>
    {
        if constexpr (new_level_information)
        {
            std::vector<coordinates_type> ret;

            index_t level_ = level();

            for (index_t i = 0; i != std::size(coordinates_); ++i)
                if (sampling_traits::is_new_point(i, level_))
                    ret.emplace_back(coordinates_[i]);

            return ret;
        }

        else
        {
            index_t level_ = level();

            auto less = Arithmetic::less(-Arithmetic::default_tolerance);
            std::set<coordinates_type, decltype(less)> set_{less};

            for (index_t level_secundum = 0; level_secundum != level_; ++level_secundum)
            {
                const sampling &sampling_ = sampling::factory(level_secundum);

                for (const coordinates_type &i : sampling_)
                    set_.insert(i);
            }

            std::vector<coordinates_type> ret;

            for (const coordinates_type &i : coordinates_)
                if (not set_.contains(i))
                    ret.emplace_back(i);

            return ret;
        }
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

