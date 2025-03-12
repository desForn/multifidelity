#pragma once

#include "arithmetic.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Print
{
    namespace Apparatus
    {
        template<class type>
        struct print_apparatus;
    }

    void print()
    {
        return;
    }

    template<class type>
    void print(const type &value);

    template<Arithmetic::index_t i = 0, class... types>
    requires(sizeof...(types) > 1)
    void print(const types &... values)
    {
        print(Utility::get_value<i>(values...));

        if constexpr (i + 1 != sizeof...(types))
            print<i + 1>(values...);

        return;
    }

    template<class... types>
    void println(const types &... values)
    {
        print(values...);
        std::cout << "\n";

        return;
    }

    namespace Apparatus
    {
        template<class type>
        struct print_apparatus
        {
            void operator()(const type &value)
            {
                std::cout << value;

                return;
            }
        };

        template<>
        struct print_apparatus<Arithmetic::complex_t>
        {
            void operator()(const Arithmetic::complex_t &value)
            {
                std::cout << "(" << value.real() << ", " << value.imag() << ")";

                return;
            }
        };

        template<class type>
        struct print_apparatus<std::vector<type>>
        {
            void operator()(const std::vector<type> &value)
            {
                std::cout << "(";

                if (!value.empty())
                    for (auto i = std::cbegin(value); i != std::cend(value); ++i)
                    {
                        print(*i);

                        if (i != std::end(value) - 1)
                            std::cout << ", ";
                    }

                std::cout << ")";

                return;
            }
        };

        template<class type, std::size_t extent>
        struct print_apparatus<std::span<type, extent>>
        {
            void operator()(const std::span<type, extent> &value)
            {
                std::cout << "(";

                if (!value.empty())
                    for (auto i = std::cbegin(value); i != std::cend(value); ++i)
                    {
                        print(*i);

                        if (i != std::end(value) - 1)
                            std::cout << ", ";
                    }

                std::cout << ")";

                return;
            }
        };

        template<class type, size_t size>
        struct print_apparatus<std::array<type, size>>
        {
            void operator()(const std::array<type, size> &value)
            {
                std::cout << "<";

                if (size != 0)
                    for (auto i = std::begin(value); i != std::end(value); ++i)
                    {
                        print(*i);

                        if (i != std::end(value) - 1)
                            std::cout << ", ";
                    }

                std::cout << ">";

                return;
            }
        };

        template<class first_type, class second_type>
        struct print_apparatus<std::pair<first_type, second_type>>
        {
            void operator()(const std::pair<first_type, second_type> &value)
            {
                std::cout << "[";
                print(std::get<0>(value));
                std::cout << ", ";
                print(std::get<1>(value));
                std::cout << "]";

                return;
            }
        };

        template<class... types>
        struct print_apparatus<std::tuple<types...>>
        {
            template<Arithmetic::index_t i = 0>
            void operator()(const std::tuple<types...> &value)
            {
                if constexpr (sizeof...(types) == 0)
                    std::cout << "[]";

                else
                {
                    if constexpr (i == 0)
                        std::cout << "[";

                    print(std::get<i>(value));

                    if constexpr (i + 1 != sizeof...(types))
                    {
                        print(", ");
                        this->operator()<i + 1>(value);
                    }

                    else
                        std::cout << "]";
                }

                return;
            }
        };
    }

    template<class type>
    void print(const type &value)
    {
        Apparatus::print_apparatus<type>{}(value);

        return;
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
