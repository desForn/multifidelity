#pragma once

#include "fwd.hpp"
#include "core/missing_data.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Invocable
{
    namespace Apparatus
    {
        class less
        {
        public:
            template<class type>
            requires(not std::floating_point<type>)
            constexpr bool operator()(const type &arg0, const type& arg1) const
                { return arg0 < arg1; }

            template<std::floating_point type>
            constexpr bool operator()(const type &arg0, const type &arg1) const
            {
                static constexpr auto l = Arithmetic::less(-Arithmetic::default_tolerance);

                return l(arg0, arg1);
            }

            template<class... types>
            constexpr bool operator()
                (const std::tuple<types...> &arg0, const std::tuple<types...> &arg1) const
            {
                return less_tuple<0, types...>(arg0, arg1);
            }

            template<index_t n, class type>
            constexpr bool operator()
                (const std::array<type, n> &arg0, const std::array<type, n> &arg1) const
            {
                for (index_t i = 0; i != n; ++i)
                {
                    if ((*this)(arg0[i], arg1[i]))
                        return true;
                    if ((*this)(arg1[i], arg0[i]))
                        return false;
                }

                return false;
            }

        private:
            template<index_t i, class... types>
            constexpr bool less_tuple
                (const std::tuple<types...> &arg0, const std::tuple<types...> &arg1) const
            {
                if ((*this)(std::get<i>(arg0), std::get<i>(arg1)))
                    return true;

                if ((*this)(std::get<i>(arg1), std::get<i>(arg0)))
                    return false;

                if constexpr (i == sizeof...(types) - 1)
                    return false;
                else
                    return less_tuple<i + 1, types...>(arg0, arg1);
            }
        };
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class output_type, class... input_types>
    class invocable_disk_binary_tree<output_type(input_types...)>
    {
        static_assert(sizeof...(input_types) != 0);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        using signature_type = output_type(input_types...);

        using key_type = std::tuple<std::remove_cvref_t<input_types>...>;
        using table_type = std::map
                <key_type, std::remove_cvref_t<output_type>, Apparatus::less>;

    private:
        using pair_type = std::pair<key_type, std::remove_cvref_t<output_type>>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        std::string data_file_;
        mutable table_type table_{};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        invocable_disk_binary_tree() = delete;
        ~invocable_disk_binary_tree() = default;

        invocable_disk_binary_tree(const invocable_disk_binary_tree &) = delete;
        invocable_disk_binary_tree &operator=(const invocable_disk_binary_tree &) = delete;

        invocable_disk_binary_tree(invocable_disk_binary_tree &&) noexcept = default;
        invocable_disk_binary_tree &operator=(invocable_disk_binary_tree &&) noexcept = default;

        invocable_disk_binary_tree(std::string data_file) :
            invocable_disk_binary_tree
                {std::move(data_file), Disk::read<pair_type>, invocable<signature_type>{}} {}
            
        template<class read_function_type>
        invocable_disk_binary_tree(std::string data_file, const read_function_type &read_function)
        requires(std::invocable<read_function_type, std::ifstream &, pair_type &>) :
            invocable_disk_binary_tree
                {std::move(data_file), read_function, invocable<signature_type>{}} {}

        invocable_disk_binary_tree(std::string data_file, invocable<signature_type>) :
            invocable_disk_binary_tree
                {std::move(data_file), Disk::read<pair_type>, invocable<signature_type>{}} {}

        template<class read_function_type>
        invocable_disk_binary_tree
            (std::string data_file, const read_function_type &read_function,
             invocable<signature_type>)
            requires(std::invocable<read_function_type, std::ifstream &, pair_type &>) :
                data_file_{std::move(data_file)}
        {        
            std::ifstream file(data_file_, std::ios::binary);
            if (not file.is_open())
                throw std::runtime_error("invocable_disk_binary_tree::constructor invalid file");

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
    public:
        template<class... input_types_bis>
        requires(is_callable_c<signature_type, input_types_bis...>)
        output_type operator()(input_types_bis &&... inputs) const
        {
            std::tuple inputs_tuple = std::make_tuple(std::forward<input_types_bis>(inputs)...);
            
            auto it = table_.find(inputs_tuple);

            if (it == std::end(table_))
                throw Core::missing_data_base{};

            return std::get<1>(*it);
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        const table_type &table() const & { return table_; }
        table_type &&table() && { return std::move(table_); }

        invocable_disk_binary_tree &clear() { table_.clear(); return *this; }
    };
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

