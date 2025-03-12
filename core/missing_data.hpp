#pragma once

#include "fwd.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Core
{
    class missing_data_base : public std::exception {};
    
    template<class container_type_>
    class missing_data : public missing_data_base
    {
    public:
        using container_type = container_type_;
        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        std::string message_{};
        container_type data_{};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        missing_data() = default;
        ~missing_data() override = default;

        missing_data(const missing_data &) = default;
        missing_data &operator=(const missing_data &) = default;

        missing_data(missing_data &&) noexcept = default;
        missing_data &operator=(missing_data &&) noexcept = default;

        missing_data(std::string message, container_type data) :
                message_{std::move(message)},
                data_{std::move(data)} {}
                
        missing_data(container_type data) :
                data_{std::move(data)} {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        [[nodiscard]] const char *what() const noexcept override { return message_.c_str(); }
        [[nodiscard]] const container_type &data() const & noexcept { return data_; }
        [[nodiscard]] container_type &data() & noexcept { return data_; }
        [[nodiscard]] container_type &&data() && noexcept { return std::move(data_); }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class function_type, class container_type>
    [[noreturn]] void throw_missing_points
        (const function_type &function, const container_type &input)
    requires(std::ranges::range<container_type> and requires
            {
                function(*(std::cbegin(input)));
            })
    {
        container_type missing_data_{};

        for (const auto &i: input)
        {
            try
            {
                (void) function(i);
            }

            catch (missing_data_base &)
            {
                if constexpr (requires(container_type c, container_type::value_type v)
                    {{c.emplace_back(v)};})
                    missing_data_.emplace_back(i);

                else if constexpr (requires(container_type c, container_type::value_type v)
                    {{c.emplace(v)};})
                    missing_data_.emplace(i);

                else
                    assert(false); // unimplemented container
            }
        }

        throw missing_data<container_type>{"Missing data", std::move(missing_data_)};
    }
}
