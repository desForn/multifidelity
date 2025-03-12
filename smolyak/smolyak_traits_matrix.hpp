#pragma once

#include "fwd.hpp"
#include "matrix/matrix_base.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Smolyak::Smolyak_traits
{
    class matrix_initialiser
    {
    public:
        using primary_domain = index_t;
        using secondary_domain = index_t;

    private:
        using matrix_type = Matrix::matrix_base<real_t>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        std::vector<std::unique_ptr<const matrix_type>> matrices_;
        std::vector<index_t> coordinates_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        matrix_initialiser() = delete;
        ~matrix_initialiser() = default;

        matrix_initialiser(const matrix_initialiser &arg) :
            coordinates_{arg.coordinates_}
        {
            matrices_.reserve(std::size(arg.matrices_));

            for (const auto &i : arg.matrices_)
                matrices_.emplace_back(i->clone());

            return;
        }

        matrix_initialiser &operator=(const matrix_initialiser &arg)
        {
            if (this == &arg)
                return *this;

            coordinates_ = arg.coordinates_;
            matrices_.clear();
            matrices_.reserve(std::size(arg.matrices_));
            for (const auto &i : arg.matrices_)
                matrices_.emplace_back(i->clone());

            return *this;
        }

        matrix_initialiser(matrix_initialiser &&) noexcept = default;
        matrix_initialiser &operator=(matrix_initialiser &&) noexcept = default;

        matrix_initialiser(const std::ranges::range auto &&matrices)
        {
            for (const matrix_type &i : matrices)
                matrices_.emplace_back(i.clone());

            ASSERT_ASSUME(not std::empty(matrices_));

            coordinates_.resize((*(std::cend(matrices_) - 1))->columns());
            std::ranges::iota(coordinates_, 0);

            ASSERT_ASSUME(std::equal(std::cbegin(matrices_) + 1, std::cend(matrices_),
                        std::begin(matrices_), [](const auto &arg0, const auto &arg1)
                        { return arg0->rows() == arg1->rows() and
                                 arg0->columns() >= arg1->columns(); }));

            ASSERT_ASSUME(std::size(coordinates_) == (*(std::cend(matrices_) - 1))->columns());

            return;
        }

        matrix_initialiser
            (const std::ranges::range auto &matrices, std::vector<index_t> coordinates) :
            coordinates_{std::move(coordinates)}
        {
            for (const matrix_type &i : matrices)
                matrices_.emplace_back(i.clone());

            ASSERT_ASSUME(not std::empty(matrices_));
            ASSERT_ASSUME(std::equal(std::cbegin(matrices_) + 1, std::cend(matrices_),
                        std::begin(matrices_), [](const auto &arg0, const auto &arg1)
                        { return arg0->rows() == arg1->rows() and
                                 arg0->columns() >= arg1->columns(); }));
            ASSERT_ASSUME(std::size(coordinates_) == (*(std::cend(matrices_) - 1))->columns());

            return;
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        std::span<const primary_domain> coordinates(index_t level) const
        {
            return {std::begin(coordinates_),
                    std::begin(coordinates_) + matrices_[level]->columns()};
        }

        std::span<const primary_domain> new_coordinates(index_t level) const
        {
            return {std::begin(coordinates_) + (level == 0 ? 0 : matrices_[level - 1]->columns()),
                    std::begin(coordinates_) + matrices_[level]->columns()};
        }

        const Matrix::matrix_base<real_t> &matrix(index_t level) const { return *matrices_[level]; }
        index_t output_size() const { return matrices_.front()->rows(); }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    
    template<bool incremental_>
    class matrix
    {
    public:
        using primary_domain = index_t;
        using secondary_domain = index_t;
        using initialiser_type = matrix_initialiser;
        static constexpr bool incremental = incremental_;

        class embedding_type;
        using dual_embedding_type = embedding_type;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        const index_t level_;
        const initialiser_type &initialiser_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        matrix() = delete;
        ~matrix() = default;

        matrix(const matrix &) = default;
        matrix &operator=(const matrix &) = default;

        matrix(matrix &&) noexcept = default;
        matrix &operator=(matrix &&) noexcept = default;

        matrix(index_t level, const initialiser_type &initialiser) :
            level_{level},
            initialiser_{initialiser} {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        void linear_operator (index_t stride, std::span<const real_t> input_span,
                std::span<real_t> output_span) const
        {
            ASSERT_ASSUME(std::size(input_span) == (std::size(coordinates()) - 1) * stride + 1);
            ASSERT_ASSUME(std::size(output_span) == (initialiser_.output_size() - 1) * stride + 1);

            const Matrix::matrix_base<real_t> &matrix = initialiser_.matrix(level_);

            for (auto it = std::cbegin(output_span); it < std::cend(output_span); it += stride)
                *it = 0;

            for (auto it = std::cbegin(matrix); it != std::cend(matrix); ++it)
                output_span[it.row() * stride] += *it * input_span[it.column() * stride];

            return;
        }
        index_t linear_operator_output_size() const { return initialiser_.output_size(); }
        real_t evaluate(std::span<const real_t> arg0, const secondary_domain &arg1) const
            { return arg0[arg1]; }
        std::span<const primary_domain> coordinates() const
            { return initialiser_.coordinates(level_); }
        std::span<const primary_domain> new_coordinates() const
            { return initialiser_.new_coordinates(level_); }
        embedding_type embedding(index_t) const;
        embedding_type dual_embedding(index_t) const;
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<bool incremental>
    class matrix<incremental>::embedding_type
    {
    private:
        const index_t output_size_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        embedding_type() = delete;
        ~embedding_type() = default;

        embedding_type(const embedding_type &) = default;
        embedding_type &operator=(const embedding_type &) = default;

        embedding_type(embedding_type &&) noexcept = default;
        embedding_type &operator=(embedding_type &&) = default;

        explicit embedding_type(index_t output_size) : output_size_{output_size} {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        void operator()(index_t, index_t, std::span<const real_t> input_span,
                std::span<real_t> output_span) const
        {
            ASSERT_ASSUME(std::size(input_span) == std::size(output_span));
            std::ranges::copy(input_span, std::begin(output_span));

            return;
        }
        index_t input_size(index_t) const { return output_size_; }
        index_t output_size(index_t) const { return output_size_; }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<bool incremental>
    matrix<incremental>::embedding_type matrix<incremental>::embedding(index_t) const
        { return embedding_type{initialiser_.output_size()}; }

    template<bool incremental>
    matrix<incremental>::embedding_type matrix<incremental>::dual_embedding(index_t) const
        { return embedding_type{initialiser_.output_size()}; }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

