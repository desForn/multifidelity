#pragma once

#include "fwd.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Invocable
{
    template<class function_type, class output_type, class... input_types>
    class invocable_function_disk_binary_tree<output_type(input_types...), function_type>
    {
        static_assert(sizeof...(input_types) != 0);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using signature_type = output_type(input_types...);

        using key_type = std::tuple<std::remove_cvref_t<input_types>...>;
        using table_type =
                std::map<key_type, std::remove_cvref_t<output_type>, std::less<key_type>>;

    private:
        using pair_type = std::pair<key_type, std::remove_cvref_t<output_type>>;

    public:
        using read_function_type = std::function<void(std::ifstream &, pair_type &)>;
        using write_function_type = std::function<void(std::ofstream &, const pair_type &)>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        const function_type function_;
        const std::string data_file_;
        const write_function_type write_function_;
        mutable table_type table_{};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        invocable_function_disk_binary_tree() = delete;
        ~invocable_function_disk_binary_tree();

        invocable_function_disk_binary_tree
                (const invocable_function_disk_binary_tree &) = delete;
        invocable_function_disk_binary_tree &operator=
                (const invocable_function_disk_binary_tree &) = delete;

        invocable_function_disk_binary_tree
                (invocable_function_disk_binary_tree &&) noexcept = default;
        invocable_function_disk_binary_tree &operator=
                (invocable_function_disk_binary_tree &&) noexcept = default;

        invocable_function_disk_binary_tree
                (function_type,
                 std::string,
                 invocable<signature_type> = {},
                 const read_function_type & = Disk::read<pair_type>,
                 write_function_type = Disk::write<pair_type>);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        template<class... input_types_bis>
        output_type operator()(input_types_bis &&...)
        const requires(is_callable_c<signature_type, input_types_bis...>);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        const function_type &function() const & { return function_; }
        function_type &&function() && { return std::move(function_); }

        const table_type &table() const & { return table_; }
        table_type &&table() && { return std::move(table_); }

        invocable_function_disk_binary_tree &clear() { table_.clear(); return *this; }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    template<class function_type, class output_type, class... input_types>
    invocable_function_disk_binary_tree<output_type(input_types...), function_type>::
    ~invocable_function_disk_binary_tree()
    {
        try
        {
            std::ofstream file(data_file_, std::ios::binary);
            if (!file.is_open())
                throw std::runtime_error
                        ("invocable_function_disk_binray_tree::destructor invalid file");

            Disk::write(file, std::size(table_));

            for (const auto &i: table_)
                write_function_(file, i);
        }

        catch (const std::exception &e)
        {
            Print::println
                    ("in invocable_function_disk_binary_tree::",
                     "~invocable_function_disk_binary_tree\n",
                     "with file name: ", data_file_,
                     ":\nException has occured when saving.\n", e.what());
        }

        return;
    }

    template<class function_type, class output_type, class... input_types>
    invocable_function_disk_binary_tree<output_type(input_types...), function_type>::
    invocable_function_disk_binary_tree
            (function_type function,
             std::string data_file,
             invocable<signature_type>,
             const read_function_type &read_function,
             write_function_type write_function) :
            function_{std::move(function)},
            data_file_{std::move(data_file)},
            write_function_{std::move(write_function)}
    {
        std::ifstream file(data_file_, std::ios::binary);

        if (!file)
            return;

        index_t n;
        Disk::read(file, n);

        pair_type node{};

        for (index_t i = 0; i != n; ++i)
        {
            read_function(file, node);
            table_.insert(node);
        }

        return;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    template<class function_type, class output_type, class... input_types>
    template<class... input_types_bis>
    output_type
    invocable_function_disk_binary_tree<output_type(input_types...), function_type>::operator()
        (input_types_bis &&...inputs) const
        requires(is_callable_c<signature_type, input_types_bis...>)
    {
        std::tuple inputs_tuple = std::make_tuple(inputs...);
        auto it = table_.find(inputs_tuple);

        if (it == std::end(table_))
        {
            pair_type node{inputs_tuple, function_(inputs...)};
            bool success;
            std::tie(it, success) = table_.insert(node);
            ASSERT_ASSUME(success);
        }

        return std::get<1>(*it);
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class function_type>
    invocable_function_disk_binary_tree(function_type, std::string) ->
    invocable_function_disk_binary_tree<invocable_signature<function_type>, function_type>;
}
