#pragma once

#include "fwd.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Invocable
{
    template<class function_type, class output_type, class input_type>
    class invocable_function_disk_sequential_vector<output_type(input_type), function_type>
    {
        static_assert(std::is_integral_v<std::remove_cvref_t<input_type>>);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using signature_type = output_type(input_type);

        using table_type = std::vector<std::remove_cvref_t<output_type>>;

        using read_function_type = std::function<void(std::ifstream &, output_type &)>;
        using write_function_type = std::function<void(std::ofstream &, const output_type &)>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        mutable table_type table_{};
        const std::string data_file_;
        function_type function_;
        const write_function_type write_function_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        invocable_function_disk_sequential_vector() = delete;

        invocable_function_disk_sequential_vector
            (const invocable_function_disk_sequential_vector &) = delete;
        invocable_function_disk_sequential_vector &operator=
            (const invocable_function_disk_sequential_vector &) = delete;

        invocable_function_disk_sequential_vector
            (invocable_function_disk_sequential_vector &&) noexcept = default;
        invocable_function_disk_sequential_vector &operator=
            (invocable_function_disk_sequential_vector &&) noexcept = default;

        invocable_function_disk_sequential_vector
                (function_type function,
                 std::string data_file,
                 read_function_type read_function = Disk::read<output_type>,
                 write_function_type write_function = Disk::write<output_type>) :
                data_file_(std::move(data_file)),
                function_(std::move(function)),
                write_function_(std::move(write_function))
        {
            std::ifstream file(data_file_, std::ios::binary);
            if (file)
            {
                index_t n;
                Disk::read(file, n);

                table_.reserve(n);

                output_type output;

                for (index_t i = 0; i != n; ++i)
                {
                    read_function(file, output);
                    table_.emplace_back(std::move(output));
                }
            }

            return;
        }

        ~invocable_function_disk_sequential_vector()
        {
            try
            {
                std::ofstream file(data_file_, std::ios::binary);

                if (!file.is_open())
                    throw std::runtime_error
                        ("invocable_function_disk_sequential_vector::destructor invalid file");

                Disk::write(file, std::size(table_));

                for (const auto &i: table_)
                    write_function_(file, i);
            }

            catch (const std::exception &e)
            {
                Print::println
                ("in invocable_function_disk_sequential_vector::",
                 "~invocable_function_disk_sequential_vector\n",
                 "with file name: ", data_file_,
                 ":\nException has occured when saving.\n", e.what());

                throw;
            }

            return;
        }

        invocable_function_disk_sequential_vector
        (const function_type &function,
         const std::string &data_file, invocable<signature_type>) :
         invocable_function_disk_sequential_vector(function, data_file) {}

        invocable_function_disk_sequential_vector(
                const function_type &function,
                const std::string &data_file,
                const read_function_type &read_function,
                const write_function_type &write_function,
                invocable<signature_type>) :
                invocable_function_disk_sequential_vector
                    (function, data_file, read_function, write_function) {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        output_type operator()(index_t input) const
        {
            if (static_cast<index_t>(input) >= std::size(table_))
            {
                table_.reserve(input + 1);

                for (index_t i = std::size(table_); i <= input; ++i)
                    table_.emplace_back(function_(i));
            }

            return table_[input];
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        const function_type &function() const & { return function_; }
        function_type &&function() && { return std::move(function_); }

        const table_type &table() const & { return table_; }
        table_type &&table() && { return std::move(table_); }

        invocable_function_disk_sequential_vector &clear() { table_.clear(); return *this; }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class function_type>
    invocable_function_disk_sequential_vector(function_type, std::string) ->
    invocable_function_disk_sequential_vector<invocable_signature<function_type>, function_type>;

    template<class function_type, class signature_type>
    invocable_function_disk_sequential_vector
        (function_type, std::string, invocable<signature_type>) ->
    invocable_function_disk_sequential_vector<signature_type, function_type>;

    template<class function_type, class read_function_type, class write_function_type>
    invocable_function_disk_sequential_vector
        (function_type, std::string, read_function_type, write_function_type) ->
    invocable_function_disk_sequential_vector<invocable_signature<function_type>, function_type>;

    template<class function_type, class read_function_type, class write_function_type, class signature_type>
    invocable_function_disk_sequential_vector
        (function_type, std::string, read_function_type, write_function_type,
         invocable<signature_type>) ->
    invocable_function_disk_sequential_vector<signature_type, function_type>;
}
