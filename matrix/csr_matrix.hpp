#pragma once

#include "matrix_base.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Matrix
{
    template<class field_type_ = real_t>
    class csr_matrix : public matrix_base<field_type_>
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
        
        struct column_value
        {
            const index_t column;
            field_type value;
        };

    private:
        using column_value_iterator_type = std::vector<column_value>::iterator;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        std::vector<column_value> column_value_{};
        std::vector<column_value_iterator_type> row_iterators_{std::end(column_value_)};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        csr_matrix() = default;
        ~csr_matrix() = default;

        csr_matrix(const csr_matrix &arg) :
            base_type{arg},
            column_value_{arg.column_value_},
            row_iterators_(rows() + 1)
        {
            std::iter_difference_t<column_value_iterator_type> offset =
                std::begin(column_value_) - std::begin(arg.column_value_);

            std::ranges::transform(arg.row_iterators_, std::begin(row_iterators_),
                    [&offset](const column_value_iterator_type &arg) -> column_value_iterator_type
                    { return arg + offset; });
            return;
        }

        csr_matrix &operator=(const csr_matrix &arg)
        {
            size() = arg.size();
            column_value_ = arg.column_value_;
            row_iterators_.clear();
            row_iterators_.resize(rows() + 1);
            
            std::iter_difference_t<column_value_iterator_type> offset =
                std::begin(column_value_) - std::begin(arg.column_value_);

            std::ranges::transform(arg.row_iterators_, std::begin(row_iterators_),
                    [&offset](const column_value_iterator_type &arg) -> column_value_iterator_type
                    { return arg + offset; });

            return *this;
        }

        csr_matrix(csr_matrix &&) noexcept = default;
        csr_matrix &operator=(csr_matrix &&) noexcept = default;

        csr_matrix(const std::array<index_t, 2> &size) : base_type{size} {}

        csr_matrix(index_t rows, index_t columns) :
            csr_matrix{std::array<index_t, 2>{rows, columns}} {}

        csr_matrix(const base_type &matrix) :
            base_type{matrix},
            row_iterators_{}
        {
            index_t previous_row = negative_1;
            auto it = std::begin(column_value_);

            for (auto i = std::cbegin(matrix); i != std::cend(matrix); ++i)
            {
                if (i.value() == 0)
                    continue;

                for(; i.row() != previous_row; ++previous_row)
                    row_iterators_.emplace_back(it);

                column_value_.emplace_back(i.column(), i.value());
                ++it;
            }
            row_iterators_.emplace_back(it);

            std::iter_difference_t<column_value_iterator_type> offset =
                std::begin(column_value_) - row_iterators_.front();

            std::ranges::transform(row_iterators_, std::begin(row_iterators_),
                    [&offset](const column_value_iterator_type &arg) -> column_value_iterator_type
                    { return arg + offset; });

            return;
        }

        csr_matrix(const std::ranges::range auto &scalars,
                const std::ranges::range auto &matrices) :
            base_type{static_cast<const base_type &>(*std::begin(matrices)).size()},
            row_iterators_{}
        {
            std::vector<matrix_iterator<const field_type>> iterators;
            std::vector<matrix_iterator<const field_type>> end_iterators;
            for (const base_type &i : matrices)
            {
                ASSERT_ASSUME(i.size() == size());
                iterators.emplace_back(std::begin(i));
                end_iterators.emplace_back(std::end(i));
            }

            index_t previous_row = negative_1;
            index_t previous_column = negative_1;
            auto it = std::begin(column_value_);

            while (true)
            {
                auto min_iterator = std::end(iterators);
                auto min_scalar = std::cend(scalars);
                
                {
                    auto i = std::begin(iterators);
                    auto j = std::begin(end_iterators);
                    auto k = std::cbegin(scalars);

                    for (; i != std::end(iterators); ++i, ++j, ++k)
                    {
                        if (*i == *j)
                            continue;

                        if (min_iterator == std::cend(iterators) or
                                (*min_iterator).row() > (*i).row() or
                                ((*min_iterator).row() == (*i).row() and
                                 (*min_iterator).column() > (*i).column()))
                        {
                            min_iterator = i;
                            min_scalar = k;
                        }
                    }
                }

                if (min_iterator == std::cend(iterators))
                    break;

                if ((*min_iterator).value() == 0)
                {
                    ++*min_iterator;
                    continue;
                }

                if (previous_row == (*min_iterator).row() and
                        previous_column == (*min_iterator).column())
                {
                    column_value_.back().value += *min_scalar * (*min_iterator).value();
                    ++*min_iterator;
                    continue;
                }

                for (; (*min_iterator).row() != previous_row; ++previous_row)
                    row_iterators_.emplace_back(it);

                column_value_.emplace_back((*min_iterator).column(), (*min_iterator).value());

                previous_column = (*min_iterator).column();
                ++it;
                ++*min_iterator;
            }
            row_iterators_.emplace_back(it);

            std::iter_difference_t<column_value_iterator_type> offset =
                std::begin(column_value_) - row_iterators_.front();

            std::ranges::transform(row_iterators_, std::begin(row_iterators_),
                    [&offset](const column_value_iterator_type &arg) -> column_value_iterator_type
                    { return arg + offset; });

            return;
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        field_type &operator[](const std::array<index_t, 2> &arg) override
        {
            ASSERT_ASSUME(arg[0] < rows() and arg[1] < columns());
            for (auto it = row_iterators_[arg[0]]; it != row_iterators_[arg[0]] + 1; ++it)
                if (it->column == arg[1])
                    return (it->value);
            throw invalid_index{arg};
        }

        iterator begin() override final;
        iterator end() override final;
        iterator row(index_t) override final;

        csr_matrix *clone() const override final { return new csr_matrix{*this}; }

        vector<field_type> operator*(const vector<field_type> &arg) const override final
        {
            ASSERT_ASSUME(std::size(arg) == columns());
            vector<field_type> ret(rows());

            auto i0 = std::begin(ret);
            auto i1 = std::begin(row_iterators_);

            for (; i0 != std::end(ret); ++i0, ++i1)
                for(auto j1 = *i1; j1 != *(i1 + 1); ++j1)
                    *i0 += j1->value * arg[j1->column];

            return ret;
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        index_t n_stored_elements() const { return std::size(column_value_); }
        index_t n_stored_elements_in_row(index_t arg) const
            { return std::distance(row_iterators_[arg + 1] - row_iterators_[arg]); }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class field_type>
        class csr_matrix_iterator_apparatus : public matrix_iterator_apparatus<field_type>
        {
        private:
            using column_value = csr_matrix<std::remove_const_t<field_type>>::column_value;
            using column_value_iterator_type = std::vector<column_value>::iterator;
            using row_iterators_bis_type = std::vector<column_value_iterator_type>::iterator;

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        private:
            column_value_iterator_type column_value_iterator_;
            row_iterators_bis_type row_iterators_begin_;
            row_iterators_bis_type row_iterators_end_;
            row_iterators_bis_type row_iterators_bis_;

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            csr_matrix_iterator_apparatus() = delete;
            ~csr_matrix_iterator_apparatus() = default;

            csr_matrix_iterator_apparatus(const csr_matrix_iterator_apparatus &) = default;
            csr_matrix_iterator_apparatus &operator=
                (const csr_matrix_iterator_apparatus &) = default;

            csr_matrix_iterator_apparatus(csr_matrix_iterator_apparatus &&) noexcept = default;
            csr_matrix_iterator_apparatus &operator=
                (csr_matrix_iterator_apparatus &&) noexcept = default;

            csr_matrix_iterator_apparatus(column_value_iterator_type column_value_iterator,
                    row_iterators_bis_type row_iterators_begin,
                    row_iterators_bis_type row_iterators_end,
                    row_iterators_bis_type row_iterators_bis) :
                column_value_iterator_{std::move(column_value_iterator)},
                row_iterators_begin_{std::move(row_iterators_begin)},
                row_iterators_end_{std::move(row_iterators_end)},
                row_iterators_bis_{std::move(row_iterators_bis)}
            {
                if (row_iterators_bis_ == row_iterators_end_)
                    return;
                while (column_value_iterator_ > *(row_iterators_bis_ + 1))
                    ++row_iterators_bis_;

                return;
            }

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            field_type &operator*() override final { return column_value_iterator_->value; }

            index_t row() const
            {
                return static_cast<index_t>
                    (std::distance(row_iterators_begin_, row_iterators_bis_));
            }

            index_t column() const { return column_value_iterator_->column; }

            csr_matrix_iterator_apparatus &operator++() override final
            {
                ++column_value_iterator_;
                if (row_iterators_bis_ == row_iterators_end_)
                    return *this;

                for (; column_value_iterator_ >= *(row_iterators_bis_ + 1); ++row_iterators_bis_)
                    ;

                return *this;
            }


            csr_matrix_iterator_apparatus operator++(int)
            {
                csr_matrix_iterator_apparatus ret{*this};
                ++*this;
                return ret;
            }

            bool operator== (const matrix_iterator_apparatus<field_type> &arg) const override final
            {
                if (typeid(*this) != typeid(arg))
                    return false;
                return column_value_iterator_ ==
                    static_cast<const csr_matrix_iterator_apparatus &>(arg).column_value_iterator_;
            }

            csr_matrix_iterator_apparatus *clone() const override final
                { return new csr_matrix_iterator_apparatus{*this}; }
        };
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class field_type>
    auto csr_matrix<field_type>::begin() -> iterator
        { return {new Apparatus::csr_matrix_iterator_apparatus<field_type>
                     (std::begin(column_value_), std::begin(row_iterators_),
                      std::end(row_iterators_), std::begin(row_iterators_))}; }

    template<class field_type>
    auto csr_matrix<field_type>::end() -> iterator
        { return {new Apparatus::csr_matrix_iterator_apparatus<field_type>
                     (std::end(column_value_), std::begin(row_iterators_),
                      std::end(row_iterators_), std::end(row_iterators_))}; }

    template<class field_type>
    auto csr_matrix<field_type>::row(index_t arg) -> iterator
        { return {new Apparatus::csr_matrix_iterator_apparatus<field_type>
                     (row_iterators_[arg], std::begin(row_iterators_),
                      std::end(row_iterators_), std::begin(row_iterators_) + arg)}; }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

