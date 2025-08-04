#pragma once

#include "vector.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Matrix
{
    template<class field_type>
    class matrix_iterator;

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class field_type_>
    class matrix_base
    {
    public:
        using field_type = field_type_;
        using iterator = matrix_iterator<field_type>;
        using const_iterator = matrix_iterator<const field_type>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        std::array<index_t, 2> size_{0, 0};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        matrix_base() = default;
        virtual ~matrix_base() = default;

        matrix_base(const matrix_base &) = default;
        matrix_base &operator=(const matrix_base &) = default;

        matrix_base(matrix_base &&) noexcept = default;
        matrix_base &operator=(matrix_base &&) noexcept = default;

        matrix_base(const std::array<index_t, 2> &arg) : size_{arg} {}
        matrix_base(index_t rows, index_t columns) : size_{rows, columns} {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        const std::array<index_t, 2> &size() const { return size_; }
        const index_t &rows() const { return size_[0]; }
        const index_t &columns() const { return size_[1]; }

        virtual field_type &operator[](const std::array<index_t, 2> &) = 0;
        const field_type &operator[](const std::array<index_t, 2> &arg) const
            { return const_cast<matrix_base &>(*this)[arg]; }
        field_type &operator[](index_t row, index_t column)
            { return (*this)[std::array<index_t, 2>{row, column}]; }
        const field_type &operator[](index_t row, index_t column) const
            { return (*this)[std::array<index_t, 2>{row, column}]; }
        
        virtual iterator begin() = 0;
        const_iterator begin() const
            { return const_iterator{std::begin(const_cast<matrix_base &>(*this))}; }
        const_iterator cbegin() const { return std::begin(*this); }

        virtual iterator end() = 0;
        const_iterator end() const
            { return const_iterator{std::end(const_cast<matrix_base &>(*this))}; }
        const_iterator cend() const { return std::end(*this); }

        virtual iterator row(index_t) = 0;
        const_iterator row(index_t arg) const
            { return const_iterator{const_cast<matrix_base &>(*this).row(arg)}; }
        const_iterator crow(index_t arg) const { return this->row(arg); }

        virtual matrix_base *clone() const = 0;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        virtual vector<field_type> operator*(const vector<field_type> &arg) const
        {
            ASSERT_ASSUME(columns() == arg.size());
            vector<field_type> ret(rows());

            for (auto i = std::cbegin(*this); i != std::cend(*this); ++i)
                ret[i.row()] += i.value() * arg[i.column()];

            return ret;
        }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    
    namespace Apparatus
    {
        template<class field_type>
        requires(not std::is_const_v<field_type>)
        class matrix_iterator_apparatus
        {
        public:
            matrix_iterator_apparatus() = default;
            virtual ~matrix_iterator_apparatus() = default;

            matrix_iterator_apparatus(const matrix_iterator_apparatus &) = default;
            matrix_iterator_apparatus &operator=(const matrix_iterator_apparatus &) = default;

            matrix_iterator_apparatus(matrix_iterator_apparatus &&) noexcept = default;
            matrix_iterator_apparatus &operator=(matrix_iterator_apparatus &&) noexcept = default;

            // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
        public:
            virtual field_type &operator*() = 0;
            virtual index_t row() const = 0;
            virtual index_t column() const = 0;
            virtual matrix_iterator_apparatus &operator++() = 0;
            virtual bool operator==(const matrix_iterator_apparatus &) const = 0;
            virtual matrix_iterator_apparatus *clone() const = 0;
        };
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class field_type_>
    class matrix_iterator
    {
        template<class field_type_secundum>
        friend class matrix_iterator;

    public:
        using field_type = field_type_;

    private:
        using non_const_field_type = std::remove_const_t<field_type>;
        using matrix_iterator_apparatus =
            Apparatus::matrix_iterator_apparatus<non_const_field_type>;
        
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        const std::unique_ptr<matrix_iterator_apparatus> pointer;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        matrix_iterator() : pointer{nullptr} {}
        ~matrix_iterator() = default;

        matrix_iterator(const matrix_iterator &arg) : pointer{arg.pointer->clone()} {}
        matrix_iterator &operator=(const matrix_iterator &arg)
        {
            pointer = arg.pointer->clone();
            return *this;
        }

        matrix_iterator(matrix_iterator &&) noexcept = default;
        matrix_iterator &operator=(matrix_iterator &&) noexcept = default;

        matrix_iterator(matrix_iterator_apparatus *arg) : pointer{arg} {}
        matrix_iterator(const matrix_iterator<std::remove_const_t<field_type>> &arg)
            requires(std::is_const_v<field_type>) : pointer{arg.pointer->clone()} {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        field_type &operator*() const { return *(*pointer); }

        field_type *operator->() { return std::addressof(*(*this)); }
        const field_type * operator->() const { return std::addressof(*(*this)); }

        index_t row() const { return pointer->row(); }
        index_t column() const { return pointer->column(); }
        decltype(auto) value(this auto &&self) { return (*self); }
        
        matrix_iterator &operator++()
        {
            ++*pointer;
            return *this;
        }

        matrix_iterator operator++(int)
        {
            matrix_iterator ret{*this};
            ++*pointer;
            return ret;
        }

        bool operator==(const matrix_iterator &arg) const { return *pointer == *arg.pointer; }
        bool operator!=(const matrix_iterator &arg) const { return not ((*this) == arg); }

        void swap(matrix_iterator &arg) { std::swap(pointer, arg.pointer); }
    };

    template<class field_type>
    void swap(matrix_iterator<field_type> &arg0, matrix_iterator<field_type> &arg1)
    {
        arg0.swap(arg1);
        return;
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

