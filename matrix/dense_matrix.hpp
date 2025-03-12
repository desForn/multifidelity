#pragma once

#include "matrix_base.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Matrix
{
    template<class field_type_ = real_t>
    class dense_matrix : public matrix_base<field_type_>
    {
    private:
        using base_type = matrix_base<field_type_>;

    public:
        using field_type = base_type::field_type;
        using iterator = base_type::iterator;
        using const_iterator = base_type::const_iterator;

        using base_type::size;
        using base_type::rows;
        using base_type::columns;
        using base_type::operator[];
        using base_type::operator*;
        using base_type::begin;
        using base_type::cbegin;
        using base_type::end;
        using base_type::cend;
        using base_type::row;
        using base_type::crow;
        
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        std::vector<field_type> values_{};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        dense_matrix() = default;
        ~dense_matrix() = default;

        dense_matrix(const dense_matrix &) = default;
        dense_matrix &operator=(const dense_matrix &) = default;

        dense_matrix(dense_matrix &&) noexcept = default;
        dense_matrix &operator=(dense_matrix &&) noexcept = default;

        dense_matrix(const std::array<index_t, 2> &size) :
            base_type{size},
            values_(size[0] * size[1], 0) {}

        dense_matrix(index_t rows, index_t columns) :
            dense_matrix{std::array<index_t, 2>{rows, columns}} {}

        dense_matrix(const std::array<index_t, 2> &size, std::vector<field_type> values) :
            base_type{size},
            values_{std::move(values)}
        { ASSERT_ASSUME(std::size(values_) == rows() * columns()); }

        dense_matrix(const base_type &matrix) :
            dense_matrix(matrix.size())
        {
            for (auto i = std::cbegin(matrix); i != std::cend(matrix); ++i)
                (*this)[i.row(), i.column()] = i.value();
            return;
        }

        dense_matrix(const std::ranges::range auto &scalars,
                const std::ranges::range auto &matrices) :
            dense_matrix(static_cast<const base_type &>(*std::begin(matrices)).size())
        {
            auto i = std::begin(scalars);
            auto j = std::begin(matrices);

            for (; i != std::end(scalars); ++i, ++j)
            {
                const auto &jm = static_cast<const base_type &>(*j);

                ASSERT_ASSUME(jm.size() == size());
                for (auto k = std::cbegin(jm); k != std::end(jm); ++k)
                    (*this)[k.row(), k.column()] += *i * k.value();
            }

            ASSERT_ASSUME(j == std::end(matrices));
            return;
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        field_type &operator[](const std::array<index_t, 2> &arg) override final
        {
            ASSERT_ASSUME(arg[0] < rows() and arg[1] < columns());
            return values_[arg[0] * columns() + arg[1] % columns()];
        }

        iterator begin() override final;
        iterator end() override final;
        iterator row(index_t) override final;

        vector<field_type> operator*(const vector<field_type> &arg) const override final
        {
            ASSERT_ASSUME(columns() == arg.size());
            vector<field_type> ret(rows());

            auto k = begin();

            for (auto &i : ret)
                for (const auto &j : arg)
                    i += *k++ * j;

            return ret;
        }

        dense_matrix *clone() const override final { return new dense_matrix{*this}; }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class field_type>
        class dense_matrix_iterator_apparatus : public matrix_iterator_apparatus<field_type>
        {
        private:
            index_t columns_;
            index_t row_{0};
            index_t column_{0};
            std::vector<field_type>::iterator value_iterator_;

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            dense_matrix_iterator_apparatus() = delete;
            ~dense_matrix_iterator_apparatus() = default;

            dense_matrix_iterator_apparatus(const dense_matrix_iterator_apparatus &) = default;
            dense_matrix_iterator_apparatus &operator=
                (const dense_matrix_iterator_apparatus &) = default;

            dense_matrix_iterator_apparatus(dense_matrix_iterator_apparatus &&) noexcept = default;
            dense_matrix_iterator_apparatus &operator=
                (dense_matrix_iterator_apparatus &&) noexcept = default;

            dense_matrix_iterator_apparatus(index_t columns,
                std::vector<field_type>::iterator value_iterator) :
                columns_{columns},
                value_iterator_{value_iterator} {}

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            field_type &operator*() override final { return *value_iterator_; }
            index_t row() const override final { return row_; }
            index_t column() const override final { return column_; }

            dense_matrix_iterator_apparatus &operator++() override final
            {
                ++column_;
                ++value_iterator_;
                if (column_ == columns_)
                {
                    column_ = 0;
                    ++row_;
                }
                return *this;
            }

            dense_matrix_iterator_apparatus operator++(int)
            {
                dense_matrix_iterator_apparatus ret{*this};
                ++*this;
                return ret;
            }

            bool operator== (const matrix_iterator_apparatus<field_type> &arg) const override final
            {
                if (typeid(*this) != typeid(arg))
                    return false;
                else
                    return value_iterator_ ==
                        static_cast<const dense_matrix_iterator_apparatus &>(arg).value_iterator_;
            }

            dense_matrix_iterator_apparatus *clone() const override final
                { return new dense_matrix_iterator_apparatus{*this}; }
        };

    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class field_type>
    auto dense_matrix<field_type>::begin() -> iterator
    {
        return {new Apparatus::dense_matrix_iterator_apparatus<field_type>
                    (columns(), std::begin(values_))};
    }

    template<class field_type>
    auto dense_matrix<field_type>::end() -> iterator
    {
        return {new Apparatus::dense_matrix_iterator_apparatus<field_type>
                    (columns(), std::end(values_))};
    }

    template<class field_type>
    auto dense_matrix<field_type>::row(index_t arg) -> iterator
    {
        return {new Apparatus::dense_matrix_iterator_apparatus<field_type>
                    (columns(), std::begin(values_) + arg * columns())};
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

