#pragma once

#include "arithmetic.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Disk
{
    namespace Apparatus
    {
        template<class type>
        struct read_apparatus;
    }

    template<class type>
    void read(std::ifstream &file, type &value)
    {
        Apparatus::read_apparatus<type>{}(file, value);

        return;
    }

    template<class first_type, class... other_types>
    requires(sizeof...(other_types) != 0)
    void read(std::ifstream &file, first_type &first_value,
        other_types &... other_values)
    {
        read(file, first_value);
        read(file, other_values...);

        return;
    } 

    namespace Apparatus
    {
        template<class type>
        struct read_apparatus
        {
            void operator()(std::ifstream &file, type &value)
            {
                file.read(reinterpret_cast<char *>(&value), sizeof(type));

                return;
            }
        };

        template<>
        struct read_apparatus<std::string>
        {
            void operator()(std::ifstream &file, std::string &value)
            {
                file >> value;

                return;
            }
        };

        template<>
        struct read_apparatus<const char *>
        {
            void operator()(std::ifstream &file, const char *&value)
            {
                ASSERT_ASSUME(value == nullptr);

                char *p = nullptr;
                char *q = nullptr;

                Arithmetic::index_t ss = 0;
                Arithmetic::index_t ps = 1;
                Arithmetic::index_t qs = 0;

                while (!file.eof())
                {
                    ps *= 2;
                    p = new char[ps];

                    ASSERT_ASSUME(p != nullptr);

                    std::memcpy(p, q, qs);
                    delete[] q;

                    while (ss < ps)
                    {
                        file.read(p + ss, 1);

                        if (p[ss] == '0')
                        {
                            value = p;
                            return;
                        }
                        ++ss;
                    }

                    q = p;

                    qs = ps;
                }

                if (ss == ps)
                {
                    q = new char[ps + 1];
                    std::memcpy(q, p, ps);
                    q[ps] = '0';
                    value = q;
                    return;
                }

                p[ps] = '0';
                value = p;
                return;
            }
        };

        template<class type>
        struct read_apparatus<std::vector<type>>
        {
            void operator()(std::ifstream &file, std::vector<type> &value)
            {
                Arithmetic::index_t size;
                read(file, size);

                value.resize(size);

                for (auto &i: value)
                    read(file, i);

                return;
            }
        };

        template<class type, size_t size>
        struct read_apparatus<std::array<type, size>>
        {
            void operator()(std::ifstream &file, std::array<type, size> &value)
            {
                for (auto &i: value)
                    read(file, i);

                return;
            }
        };

        template<class first_type, class second_type>
        struct read_apparatus<std::pair<first_type, second_type>>
        {
            void operator()(std::ifstream &file, std::pair<first_type, second_type> &value)
            {
                read(file, std::get<0>(value));
                read(file, std::get<1>(value));

                return;
            }
        };

        template<class... types>
        struct read_apparatus<std::tuple<types...>>
        {
            template<Arithmetic::index_t i = 0>
            void operator()(std::ifstream &file, std::tuple<types...> &value)
            {
                static_assert(sizeof...(types) != 0);

                read(file, std::get<i>(value));

                if constexpr (i + 1 != sizeof...(types))
                    this->operator()<i + 1>(file, value);

                return;
            }
        };
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    namespace Apparatus
    {
        template<class type>
        struct write_apparatus;
    }

    template<class type>
    void write(std::ofstream &file, const type &value)
    {
        Apparatus::write_apparatus<type>{}(file, value);

        return;
    }

    template<class first_type, class... other_types>
    void write(std::ofstream &file, const first_type &first_value,
        const other_types &... other_values)
    {
        write(file, first_value);
        write(file, other_values...);

        return;
    } 

    namespace Apparatus
    {
        template<class type>
        struct write_apparatus
        {
            void operator()(std::ofstream &file, const type &value)
            {
                file.write(reinterpret_cast<const char *>(&value), sizeof(type));

                return;
            }
        };

        template<>
        struct write_apparatus<std::string>
        {
            void operator()(std::ofstream &file, const std::string &value)
            {
                file << value;

                return;
            }
        };

        template<>
        struct write_apparatus<const char *>
        {
            void operator()(std::ofstream &file, const char * const &value)
            {
                file.write(value, std::strlen(value));

                return;
            }
        };

        template<std::size_t n>
        struct write_apparatus<const char[n]>
        {
            void operator()(std::ofstream &file, const char (&value) [n])
            {
                write_apparatus<const char *>(file, value);
                return;
            }
        };

        template<>
        struct write_apparatus<char *>
        {
            void operator()(std::ofstream &file, char * const &value)
            {
                write_apparatus<const char *>{}(file, value);
                return;
            }
        };

        template<std::size_t n>
        struct write_apparatus<char[n]>
        {
            void operator()(std::ofstream &file, const char (&value) [n])
            {
                write_apparatus<const char *>{}(file, value);
                return;
            }
        };

        template<class type>
        struct write_apparatus<std::vector<type>>
        {
            void operator()(std::ofstream &file, const std::vector<type> &value)
            {
                Arithmetic::index_t size = std::size(value);
                write(file, size);

                for (const auto &i: value)
                    write(file, i);

                return;
            }
        };

        template<class type, size_t size>
        struct write_apparatus<std::array<type, size>>
        {
            void operator()(std::ofstream &file, const std::array<type, size> &value)
            {
                for (const auto &i: value)
                    write(file, i);

                return;
            }
        };

        template<class first_type, class second_type>
        struct write_apparatus<std::pair<first_type, second_type>>
        {
            void operator()(std::ofstream &file, const std::pair<first_type, second_type> &value)
            {
                write(file, std::get<0>(value));
                write(file, std::get<1>(value));

                return;
            }
        };

        template<class... types>
        struct write_apparatus<std::tuple<types...>>
        {
            template<Arithmetic::index_t i = 0>
            void operator()(std::ofstream &file, const std::tuple<types...> &value)
            {
                static_assert(sizeof...(types) != 0);

                write(file, std::get<i>(value));

                if constexpr (i + 1 != sizeof...(types))
                    this->operator()<i + 1>(file, value);

                return;
            }
        };
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class type, std::size_t n>
    void write_array(const std::string &file_name, const std::array<type, n> &array)
    {
        try
        {
            std::ofstream data_file(file_name, std::ios::binary);
            for (std::size_t i = 0; i != n; ++i)
                data_file.write(reinterpret_cast<const char *>(&array[i]), sizeof(type));

            data_file.close();
        }

        catch (std::exception &e)
        {
            throw;
        }

        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class type, std::size_t n>
    void write_array(std::ofstream &data_file, const std::array<type, n> &array)
    {
        for (std::size_t i = 0; i != n; ++i)
            data_file.write(reinterpret_cast<const char *>(&array[i]), sizeof(type));

        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class type, std::size_t n>
    std::array<type, n> write_array(std::string &file_name)
    {
        std::ofstream data_file(file_name, std::ios::binary);

        ASSERT_ASSUME(data_file);

        std::array<type, n> array;
        write_array(data_file, array);

        data_file.close();

        return array;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class type, std::size_t n>
    void read_array(const std::string &file_name, std::array<type, n> &array)
    {
        try
        {
            std::ifstream data_file(file_name, std::ios::binary);
            for (std::size_t i = 0; i != n; ++i)
                data_file.read(reinterpret_cast<char *>(&array[i]), sizeof(type));

            data_file.close();
        }

        catch (std::exception &e)
        {
            throw;
        }

        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class type, std::size_t n>
    std::array<type, n> read_array(const std::string &file_name)
    {
        std::array<type, n> array;
        read_array(file_name, array);

        return array;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class type, std::size_t n>
    void read_array(std::ifstream &data_file, std::array<type, n> &array)
    {
        for (std::size_t i = 0; i != n; ++i)
            data_file.read(reinterpret_cast<char *>(&array[i]), sizeof(type));

        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class type, std::size_t n>
    std::array<type, n> read_array(std::ifstream &data_file)
    {
        std::array<type, n> array;
        read_array(data_file, array);

        return array;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class type, std::size_t n>
    std::array<type, n> read_array(std::string &file_name)
    {
        std::ifstream data_file(file_name, std::ios::binary);

        ASSERT_ASSUME(data_file);

        std::array<type, n> array;
        read_array(data_file, array);

        data_file.close();

        return array;
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
