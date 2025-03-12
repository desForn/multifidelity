#pragma once

#include "fwd.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Matrix
{
    template<class field_type_ = real_t>
    class vector
    {
    public:
        using field_type = field_type_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        std::vector<field_type> vector_{};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        vector() = default;
        ~vector() = default;

        vector(const vector &) = default;
        vector &operator=(const vector &) = default;

        vector(vector &&) noexcept = default;
        vector &operator=(vector &&) noexcept = default;

        vector(index_t size) : vector_(size, 0) {}
        
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        index_t size() const { return std::size(vector_); }
        decltype(auto) operator[](this auto &&self, index_t arg) { return (self.vector_[arg]); }

        auto begin(this auto &&self) { return std::begin(self.vector_); }
        auto cbegin(this auto &&self) { return std::cbegin(self.vector_); }
        auto end(this auto &&self) { return std::end(self.vector_); }
        auto cend(this auto &&self) { return std::cend(self.vector_); }
    };
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

