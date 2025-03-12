#pragma once

#include "fwd.hpp"
#include "invocable/evaluate.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Smolyak
{
    enum verbosity : index_t
    {
        levels_prior_adaptivity = 1,
        forward_neighborhood = 2,
        empty_forward_neighborhood = 4,
        information_of_each_level = 8,
        new_level = 16,
        none = 0,
        all = negative_1,
    };

    enum class node_state
    {
        unfitted,       // inactive and not fitted
        active,         //   active and     fitted
        inactive        // inactive and     fitted
    };

    template<class function_t, class cost_function_t, smolyak_traits_c... traits_types>
    class smolyak_approximation
    {
    public:
        using function_type = function_t;
        using cost_function_type = cost_function_t;
        static constexpr index_t n_variates = sizeof...(traits_types);

        using traits_type = std::tuple<traits_types...>;
        using level_type = std::array<index_t, n_variates>;

    private:
        static constexpr std::array<index_t, n_variates> incremental_
                {(traits_types::incremental ? 1 : 2)...};

        using secondary_coordinate_type_apparatus =
            Utility::non_void_tuple<typename traits_types::secondary_domain...>;
        using traits_initialiser_type_apparatus =
            Utility::non_void_tuple<typename traits_types::initialiser_type...>;

    public:
        class void_type {};

        using primary_coordinate_type = std::tuple<typename traits_types::primary_domain...>;

        static constexpr bool void_secondary_domain =
            std::same_as<secondary_coordinate_type_apparatus, void>;
        using secondary_coordinate_type =
                std::conditional_t<void_secondary_domain, void_type,
                                   secondary_coordinate_type_apparatus>;

        static constexpr bool void_initialiser =
            std::same_as<traits_initialiser_type_apparatus, void>;
        using traits_initialiser_type =
            std::conditional_t<void_initialiser, void_type, traits_initialiser_type_apparatus>;

        static constexpr bool defined_inner_product =
            (Apparatus::smolyak_traits_inner_product_c<traits_types> and ...);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        struct adaptivity_information
        {
            level_type added_level;
            real_t node_norm;
            real_t error_estimation;
        };

    private:
        class operator_type;
        using level_set_type = std::unordered_set<level_type, Core::hash>;
        using operator_map_type = std::unordered_map<level_type, operator_type, Core::hash>;
        using inner_product_map_type =
            std::unordered_map<std::array<level_type, 2>, real_t, Core::hash>;

    public:
        using coordinates_set_type = std::unordered_set<primary_coordinate_type, Core::hash>;
        using missing_data = Core::missing_data<coordinates_set_type>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        const function_type function_;
        const cost_function_type cost_function_;
        const traits_initialiser_type traits_initialiser_;
        const level_type traits_maximum_level_;
        level_type maximum_level_{Utility::uniform_array<index_t, n_variates>(negative_1)};
        operator_map_type nodes_{};
        mutable inner_product_map_type inner_product_{};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        smolyak_approximation() = delete;
        ~smolyak_approximation() = default;

        smolyak_approximation(const smolyak_approximation &) = delete;
        smolyak_approximation &operator=(const smolyak_approximation &) = delete;

        smolyak_approximation(smolyak_approximation &&) noexcept = default;
        smolyak_approximation &operator=(smolyak_approximation &&) noexcept = default;

        explicit smolyak_approximation
            (function_type arg0, cost_function_type arg1) requires(void_initialiser) :
                function_{std::move(arg0)},
                cost_function_{std::move(arg1)},
                traits_initialiser_{},
                traits_maximum_level_{traits_maximum_level<n_variates - 1>()} {}

        explicit smolyak_approximation
            (function_type arg0, cost_function_type arg1, traits_initialiser_type arg2)
            requires(not void_initialiser) :
                function_{std::move(arg0)},
                cost_function_{std::move(arg1)},
                traits_initialiser_{std::move(arg2)},
                traits_maximum_level_{traits_maximum_level
                    <n_variates - 1, std::tuple_size_v<traits_initialiser_type> - 1>
                        (traits_initialiser_)} {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        template<index_t i>
        static std::array<index_t, i + 1> traits_maximum_level() requires(void_initialiser);

        template<index_t i, index_t j>
        static std::array<index_t, i + 1> traits_maximum_level
            (const traits_initialiser_type &) requires(not void_initialiser);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        [[nodiscard]] const function_type &function() const { return function_; }
        [[nodiscard]] const cost_function_type &cost_function() const { return cost_function_; }
        [[nodiscard]] decltype(auto) maximum_level(this auto&& self)
            { return (self.maximum_level_); }
        [[nodiscard]] const operator_map_type &nodes() const { return nodes_; }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        real_t operator()() const requires(void_secondary_domain);
        real_t operator()(const secondary_coordinate_type &) const
            requires(not void_secondary_domain);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        [[nodiscard]] std::set<primary_coordinate_type, decltype(Arithmetic::less())>
            coordinates() const;
        [[nodiscard]] real_t norm() const;
        [[nodiscard]] real_t cost() const;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        [[nodiscard]] real_t inner_product(const level_type &, const level_type &) const;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        void fit_node_apparatus(const level_type &);
        void activate_node_apparatus(const level_type &);
        void deactivate_node_apparatus(const level_type &);
        bool valid_to_spawn(const level_type &arg) const
        {
            return std::ranges::equal(arg, maximum_level_, std::less_equal{}) and
                std::ranges::equal(arg, traits_maximum_level_, std::less_equal{});
        }

    public:
        smolyak_approximation &fit_node(const level_type &arg)
            { return fit_nodes(level_set_type{arg}); }
        smolyak_approximation &fit_nodes(const Traits::range_c<level_type> auto &arg)
            { return fit_nodes(level_set_type{arg}); }
        smolyak_approximation &fit_nodes(level_set_type);
        smolyak_approximation &activate_node(const level_type &arg)
            { return activate_nodes(level_set_type{arg}); }
        smolyak_approximation &activate_nodes(const Traits::range_c<level_type> auto &arg)
            { return activate_nodes(level_set_type{arg}); }
        smolyak_approximation &activate_nodes(level_set_type);
        smolyak_approximation &deactivate_node(const level_type &arg)
            { return deactivate_nodes(level_set_type{arg}); }
        smolyak_approximation &deactivate_nodes(const Traits::range_c<level_type> auto &arg)
            { return deactivate_nodes(level_set_type{arg}); }
        smolyak_approximation &deactivate_nodes(level_set_type);

        smolyak_approximation &activate_nodes_by_total_order(index_t);
        adaptivity_information activate_node_adaptively(index_t = 0)
        requires(defined_inner_product);

        smolyak_approximation &activate_all_nodes();
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<smolyak_traits_c... traits_types>
    struct smolyak_approximation_factory
    {
        template<class function_type, class cost_function_type>
        [[nodiscard]] static smolyak_approximation
                <std::remove_reference_t<function_type>,
                 std::remove_reference_t<cost_function_type>, traits_types...>
        factory(function_type &&function, cost_function_type &&cost_function)
        {
            return smolyak_approximation<std::remove_reference_t<function_type>,
                                         std::remove_reference_t<cost_function_type>,
                                         traits_types...>
                    {std::forward<function_type>(function),
                     std::forward<cost_function_type>(cost_function)};
        }

        template<class function_type, class cost_function_type,
                 class... constructor_arguments_types>
        [[nodiscard]] static smolyak_approximation
                <std::remove_reference_t<function_type>,
                 std::remove_reference_t<cost_function_type>, traits_types...>
        factory(function_type &&function, cost_function_type &&cost_function,
                constructor_arguments_types &&... constructor_arguments)
        {
            return smolyak_approximation<std::remove_reference_t<function_type>,
                                         std::remove_reference_t<cost_function_type>,
                                         traits_types...>
                    {std::forward<function_type>(function),
                     std::forward<cost_function_type>(cost_function),
                     std::forward<constructor_arguments_types>(constructor_arguments)...};
        }
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    class smolyak_approximation<function_type, cost_function_type, traits_types...>::operator_type
    {
    private:
        const level_type level_;
        const traits_type traits_;
        integer_t coefficient_{0};
        real_t cost_{-1};
        real_t norm_{-1};
        node_state node_state_{node_state::unfitted};
        std::vector<real_t> v_{};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        operator_type() = delete;
        ~operator_type() = default;

        operator_type(const operator_type &) = delete;
        operator_type &operator=(const operator_type &) = delete;

        operator_type(operator_type &&) noexcept = default;
        operator_type &operator=(operator_type &&) noexcept = default;

        explicit operator_type(const level_type &arg0) requires(void_initialiser) :
                level_{arg0},
                traits_{level_} {}
        explicit operator_type(const level_type &arg0, const traits_initialiser_type &arg1)
            requires(not void_initialiser) :
                level_{arg0},
                traits_{traits_factory(level_, arg1)} {}

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        template<index_t i = n_variates - 1,
                 index_t j = std::tuple_size_v<traits_initialiser_type> - 1>
        static auto traits_factory (const level_type &, const traits_initialiser_type &)
            requires(not void_initialiser);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        const level_type &level() const { return level_; }
        decltype(auto) coefficient(this auto &&self) { return (self.coefficient_); }
        decltype(auto) state(this auto &&self) { return (self.node_state_); }
        const real_t &cost() const
        {
            ASSERT_ASSUME(node_state_ != node_state::unfitted);
            ASSERT_ASSUME(cost_ != -1);
            return cost_;
        }
        decltype(auto) norm(this auto &&self)
        {
            ASSERT_ASSUME(self.node_state_ != node_state::unfitted);
            return (self.norm_);
        }
        const traits_type &traits() const { return traits_; }
        std::span<const real_t> v() const
        {
            ASSERT_ASSUME(node_state_ != node_state::unfitted);
            return v_;
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        operator_type &fit(const function_type &, const cost_function_type &);
        [[nodiscard]] real_t evaluate() const requires(void_secondary_domain);
        [[nodiscard]] real_t evaluate(const secondary_coordinate_type &) const
        requires(not void_secondary_domain);
        [[nodiscard]] std::vector<real_t> homomorphism(const homomorphism_c auto &...homomorphisms)
        const requires(sizeof...(homomorphisms) == n_variates);
        [[nodiscard]] real_t inner_product(const operator_type &) const;

    private:
        template<index_t i_variate = n_variates - 1>
        void apply_linear_operator(std::vector<real_t> &, std::vector<real_t> &, index_t);

        template<index_t i_variate = n_variates - 1>
        void evaluate(std::vector<real_t> &, std::vector<real_t> &) const
        requires(void_secondary_domain);

        template<index_t i_variate = n_variates - 1,
                index_t j_variate = std::tuple_size_v<secondary_coordinate_type> - 1>
        void evaluate(std::vector<real_t> &, std::vector<real_t> &,
                const secondary_coordinate_type &) const requires(not void_secondary_domain);

        template<index_t i_variate = n_variates - 1>
        void homomorphism(std::vector<real_t> &, std::vector<real_t> &, index_t,
                 const homomorphism_c auto &...homomorphisms) const
        requires(sizeof...(homomorphisms) == n_variates);

        template<index_t i_variate = n_variates - 1>
        index_t homomorphism_apparatus_maximum_size
            (index_t, index_t, const homomorphism_c auto &... homomorphisms) const
        requires(sizeof...(homomorphisms) == n_variates);
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    template<index_t i>
    std::array<index_t, i + 1> smolyak_approximation
        <function_type, cost_function_type, traits_types...>::traits_maximum_level()
        requires(void_initialiser)
    {
        static_assert(i < n_variates);

        using traits_type = Utility::get_type<i, traits_types...>;
        using ret_type = std::array<index_t, i + 1>;

        constexpr bool max_level_c = requires(const Utility::get_type<i, traits_types...> traits)
        {{traits.max_level()}; };

        if constexpr (i == 0)
        {
            if constexpr (max_level_c)
                return ret_type{traits_type{0}.max_level()};
            else
                return ret_type{negative_1};
        }
        else
        {
            std::array<index_t, i> recursive_ret = traits_maximum_level<i - 1>();
            ret_type ret;
            std::copy(std::cbegin(recursive_ret), std::cend(recursive_ret), std::begin(ret));

            if constexpr (max_level_c)
                ret.back() = traits_type{0}.max_level();
            else
                ret.back() = negative_1;
            return ret;
        }
    }
    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    template<index_t i, index_t j>
    std::array<index_t, i + 1> smolyak_approximation
        <function_type, cost_function_type, traits_types...>::traits_maximum_level
        (const traits_initialiser_type &arg) requires(not void_initialiser)
    {
        static_assert(i < n_variates and j < std::tuple_size_v<traits_initialiser_type> and
            (i != 0 or j == 0));

        using traits_type = Utility::get_type<i, traits_types...>;
        using ret_type = std::array<index_t, i + 1>;

        constexpr bool max_level_c = requires(const Utility::get_type<i, traits_types...> traits)
        {{traits.max_level()}; };

        if constexpr (std::same_as<typename traits_type::initialiser_type, void>)
        {
            if constexpr (max_level_c)
            {
                traits_type traits{0};

                if constexpr (i == 0)
                    return ret_type{traits.max_level()};
                else
                {
                    std::array<index_t, i> recursive_ret = traits_maximum_level<i - 1, j>(arg);
                    ret_type ret;
                    std::copy(std::cbegin(recursive_ret), std::cend(recursive_ret),
                              std::begin(ret));

                    ret.back() = traits.max_level();
                    return ret;
                }
            }
            else
            {
                if constexpr (i == 0)
                    return ret_type{negative_1};
                else
                {
                    std::array<index_t, i> recursive_ret = traits_maximum_level<i - 1, j>(arg);
                    ret_type ret;
                    std::copy(std::cbegin(recursive_ret), std::cend(recursive_ret),
                              std::begin(ret));

                    ret.back() = negative_1;
                    return ret;
                }
            }
        }
        else
        {
            if constexpr (max_level_c)
            {
                traits_type traits{0, std::get<j>(arg)};

                if constexpr (i == 0)
                    return ret_type{traits.max_level()};
                else
                {
                    std::array<index_t, i> recursive_ret = traits_maximum_level<i - 1, j - 1>(arg);
                    ret_type ret;
                    std::copy(std::cbegin(recursive_ret), std::cend(recursive_ret),
                              std::begin(ret));

                    ret.back() = traits.max_level();
                    return ret;
                }
            }
            else
            {
                if constexpr (i == 0)
                    return ret_type{negative_1};
                else
                {
                    std::array<index_t, i> recursive_ret = traits_maximum_level<i - 1, j - 1>(arg);
                    ret_type ret;
                    std::copy(std::cbegin(recursive_ret), std::cend(recursive_ret),
                              std::begin(ret));

                    ret.back() = negative_1;
                    return ret;
                }
            }
        }
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    real_t smolyak_approximation<function_type, cost_function_type, traits_types...>::operator()()
        const requires(void_secondary_domain)
    {
        real_t ret{0};

        for (const auto &it: nodes_)
            if (it.second.coefficient())
                ret += it.second.coefficient() * it.second.evaluate();

        return ret;
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    real_t smolyak_approximation<function_type, cost_function_type, traits_types...>::operator()
        (const secondary_coordinate_type &arg) const requires(not void_secondary_domain)
    {
        real_t ret{0};

        for (const auto &it: nodes_)
            if (it.second.coefficient())
                ret += it.second.coefficient() * it.second.evaluate(arg);

        return ret;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    auto smolyak_approximation<function_type, cost_function_type, traits_types...>::
        coordinates() const -> std::set<primary_coordinate_type, decltype(Arithmetic::less())>
    {
        std::set<primary_coordinate_type, decltype(Arithmetic::less())> ret;

        auto carteisan_product = []<class... traits_types_secundum>
            (const traits_types_secundum &... traits)
        {
            return Utility::cartesian_product(traits.coordinates()...);
        };

        for (const auto &i: nodes_)
        {
            std::vector<primary_coordinate_type> c = std::apply(carteisan_product, i.traits());
            ret.insert_range(c);
        }
        
        return ret;
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    real_t smolyak_approximation<function_type, cost_function_type, traits_types...>::norm() const
    {
        real_t ret{0};
        
        for (const auto &i : nodes_)
        {
            if (i.second.state() == node_state::unfitted)
                continue;

            for (const auto &j : nodes_)
            {
                if (j.second.state() == node_state::unfitted)
                    continue;

                if (j.first < i.first)
                    continue;

                real_t a = inner_product(i.first, j.first) *
                        i.second.coefficient() * j.second.coefficient();

                if (j.first != i.first)
                    a *= 2;
                ret += a;
            }
        }

        ASSERT_ASSUME(Arithmetic::less()(0, ret));
        if (ret > 0) [[likely]]
            return std::sqrt(ret);
        return 0;
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    real_t smolyak_approximation<function_type, cost_function_type, traits_types...>::cost() const
    {
        real_t ret{0};

        for (const auto &node : nodes_)
            if (node.second.state() != node_state::unfitted)
                ret += node.second.cost();

        return ret;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    real_t smolyak_approximation<function_type, cost_function_type, traits_types...>::
        inner_product(const level_type &arg0, const level_type &arg1) const
    {
        ASSERT_ASSUME(arg0 <= arg1);

        std::array<level_type, 2> key{arg0, arg1};

        auto it = inner_product_.find(key);
        if (it != std::end(inner_product_))
            return it->second;

        auto it0 = nodes_.find(arg0);
        auto it1 = nodes_.find(arg1);

        ASSERT_ASSUME(it0 != std::end(nodes_) and it1 != std::end(nodes_));

        real_t ret = it0->second.inner_product(it1->second);

        inner_product_[key] = ret;
        return ret;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    void smolyak_approximation<function_type, cost_function_type, traits_types...>::
        fit_node_apparatus(const level_type &arg)
    {
        auto node = nodes_.find(arg);
        if (node == std::end(nodes_))
        {
            if constexpr (void_initialiser)
                node = nodes_.emplace(arg, operator_type{arg}).first;
            else
                node = nodes_.emplace(arg, operator_type{arg, traits_initialiser_}).first;
        }

        if (node->second.state() == node_state::unfitted)
            node->second.fit(function_, cost_function_); // May throw missing_data

        std::unordered_map<level_type, integer_t, Core::hash> coefficients_predecessors;
        for (const level_type &a : Utility::counter{incremental_})
        {
            if (not std::ranges::equal
                    (arg, a, [](index_t arg0, index_t arg1) { return arg0 >= arg1; }))
                continue;

            level_type al = Utility::subtract_array(arg, a);
            index_t parity = std::accumulate(std::cbegin(a), std::cend(a), 0) % 2;
            integer_t c = parity ? -1 : 1;
            coefficients_predecessors[al] = c;
        }

        real_t norm = 0;
        for (const auto &i : coefficients_predecessors)
            for (const auto &j : coefficients_predecessors)
            {
                if (i.first > j.first)
                    continue;
                real_t c = inner_product(i.first, j.first) * i.second * j.second;
                if (i.first == j.first)
                    norm += c;
                else
                    norm += 2 * c;
            }

        ASSERT_ASSUME(Arithmetic::less()(0, norm));
        if (norm < 0)
            node->second.norm() = 0;
        node->second.norm() = std::sqrt(norm);

        node->second.state() = node_state::inactive;

        return;
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    void smolyak_approximation<function_type, cost_function_type, traits_types...>::
        activate_node_apparatus(const level_type &arg)
    {
        operator_type &node = nodes_.find(arg)->second;

        if (node.state() == node_state::active)
            return;

        for (const level_type &a : Utility::counter(incremental_))
        {
            if (not std::ranges::equal
                    (arg, a, [](index_t arg0, index_t arg1) { return arg0 >= arg1; }))
                continue;

            index_t parity = std::accumulate(std::cbegin(a), std::cend(a), 0) % 2;
            nodes_.find(Utility::subtract_array(arg, a))->second.coefficient() += parity ? -1 : 1;
        }

        node.state() = node_state::active;

        return;
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    void smolyak_approximation<function_type, cost_function_type, traits_types...>::
        deactivate_node_apparatus(const level_type &arg)
    {
        operator_type &node = nodes_.find(arg).second;

        if (node.state() != node_state::active)
            return;

        for (const level_type &a : Utility::counter(incremental_))
        {
            if (not std::ranges::equal
                    (arg, a, [](index_t arg0, index_t arg1) { return arg0 >= arg1; }))
                continue;

            index_t parity = std::accumulate(std::begin(a), std::end(a), 0) % 2;
            nodes_.find(Utility::subtract_array(arg, a))->second.coefficient() += parity ? 1 : -1;
        }

        node.state() = node_state::inactive;

        return;
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    auto smolyak_approximation<function_type, cost_function_type, traits_types...>::fit_nodes
        (level_set_type arg) -> smolyak_approximation &
    {
        /* If no exception is thown after execution of the function all nodes in arg and their
         * predecessors are fitted. The state of the rest of nodes in not changed.
         * This function may throw missing_data
         *      It has a particular strong exception guarantee.
         *      If thrown, the node state of all nodes is the same as before calling. However, all
         *      requested levels and its predecessors will be equally added but set to unfitted or
         *      inactive. */

        for (auto it = std::begin(arg); it != std::end(arg);)
        {
            if (not valid_to_spawn(*it))
            {
                it = arg.erase(it);
                continue;
            }

            auto node_it = nodes_.find(*it);
            if (node_it != std::end(nodes_) and
                node_it->second.state() != node_state::unfitted)
                it = arg.erase(it);
            else
                ++it;
        }

        level_set_type unfitted_predecessors;

        for (const level_type &level : arg)
            for (index_t i = 0; i != n_variates; ++i)
            {
                if (level[i] == 0)
                    continue;

                level_type predecessor = level;
                --predecessor[i];

                auto node_it = nodes_.find(predecessor);
                if (node_it == std::end(nodes_) or node_it->second.state() == node_state::unfitted)
                    unfitted_predecessors.emplace(predecessor);
            }

        missing_data missing_data_;

        if (not std::empty(unfitted_predecessors))
        {
            try
            {
                fit_nodes(std::move(unfitted_predecessors));
            }

            catch (missing_data &e)
            {
                missing_data_.data().merge(e.data());
            }
        }

        for (const level_type &level : arg)
        {
            try
            {
                fit_node_apparatus(level);
            }

            catch (missing_data &e)
            {
                missing_data_.data().merge(e.data());
            }
        }

        if (not std::empty(missing_data_.data()))
            throw missing_data_;

        return *this;
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    auto smolyak_approximation<function_type, cost_function_type, traits_types...>::activate_nodes
        (level_set_type arg) -> smolyak_approximation &
    {
        fit_nodes(arg);

        for (auto it = std::begin(arg); it != std::end(arg);)
        {
            if (not valid_to_spawn(*it))
            {
                it = arg.erase(it);
                continue;
            }

            if (nodes_.find(*it)->second.state() == node_state::active)
                it = arg.erase(it);
            else
                ++it;
        }

        level_set_type inactive_predecessors;

        for (const level_type &level : arg)
            for (index_t i = 0; i != n_variates; ++i)
            {
                if (level[i] == 0)
                    continue;

                level_type predecessor = level;
                --predecessor[i];

                if (nodes_.find(predecessor)->second.state() == node_state::inactive)
                    inactive_predecessors.emplace(predecessor);
            }

        if (not std::empty(inactive_predecessors))
            activate_nodes(std::move(inactive_predecessors));

        for (const level_type &level : arg)
            activate_node_apparatus(level);

        return *this;
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    auto smolyak_approximation<function_type, cost_function_type, traits_types...>::deactivate_nodes
            (level_set_type arg) -> smolyak_approximation &
    {
        for (auto it = std::begin(arg); it != std::end(arg);)
        {
            auto node_it = nodes_.find(*it);
            if (node_it == std::end(nodes_) or node_it->second.state() != node_state::active)
                it = arg.erase(it);
            else
                ++it;
        }

        level_set_type active_successors;

        for (const level_type &level : arg)
            for (index_t i = 0; i != n_variates; ++i)
            {
                level_type successor = level;
                ++successor[i];

                auto node_it = nodes_.find(successor);
                if (node_it != std::end(nodes_) and
                    node_it->second.state() == node_state::active)
                    active_successors.emplace(successor);
            }

        if (not std::empty(active_successors))
            deactivate_nodes(std::move(active_successors));

        for (const level_type &level : arg)
            deactivate_node_apparatus(level);

        return *this;
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    auto smolyak_approximation<function_type, cost_function_type, traits_types...>::
        activate_nodes_by_total_order(index_t arg) -> smolyak_approximation &
    {
        level_set_type levels;

        level_type i{};
        index_t order = 0;

        while (true)
        {
            levels.emplace(i);

            for (index_t j = n_variates - 1; j != negative_1; --j)
            {
                if (order != arg)
                {
                    ++order;
                    ++i[j];
                    goto add_nodes_by_total_order_next_iteration;
                }
                order -= i[j];
                i[j] = 0;
            }
            break;
            add_nodes_by_total_order_next_iteration:
        }

        activate_nodes(std::move(levels));
        return *this;
    };

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    auto smolyak_approximation<function_type, cost_function_type, traits_types...>::
        activate_node_adaptively(index_t arg) -> adaptivity_information
    requires(defined_inner_product)
    {
        static constexpr adaptivity_information error_ret
            {Utility::uniform_array<index_t, n_variates>(negative_1), -1, -1};

        if (std::empty(nodes_))
        {
            level_type level{};
            if (arg & verbosity::levels_prior_adaptivity)
                Print::println("No levels prior adaptivity.");

            if (not valid_to_spawn(level))
            {
                if (arg & verbosity::empty_forward_neighborhood)
                    Print::println("Empty forward neighborhood");

                return error_ret;
            }
            if (arg & verbosity::forward_neighborhood)
            {
                Print::println("Forward_neighborhood:");
                Print::println("\t", level);
            }

            fit_node(level);

            auto it = nodes_.find(level);

            real_t norm = it->second.norm();
            real_t cost = it->second.cost();

            if (arg & verbosity::information_of_each_level)
            {
                Print::println("Information of each level:");
                Print::println("\tLevel = ", level, "; Error = ", norm, "; Cost = ", cost,
                               "; ratio = ", norm / cost);
            }

            if (arg & verbosity::new_level)
                Print::println("New level = ", level);

            activate_node(level);
            return {level, norm, norm};
        }

        if (arg & verbosity::levels_prior_adaptivity)
        {
            Print::println("Levels prior adaptivity:");
            for (const auto &node: nodes_)
                if (node.second.state() == node_state::active)
                    Print::println("\t", node.first);
        }

        level_set_type forward_neighborhood;
        for (const auto &node : nodes_)
        {
            if (node.second.state() != node_state::active)
                continue;

            for (index_t i = 0; i != n_variates; ++i)
            {
                level_type a = Utility::increment_copy(node.first, i);

                if (not valid_to_spawn(a))
                    continue;

                auto node_it = nodes_.find(a);
                if (node_it != std::end(nodes_) and
                    node_it->second.state() == node_state::active)
                    continue;

                for (index_t j = 0; j != n_variates; ++j)
                {
                    if (i == j or a[j] == 0)
                        continue;

                    level_type b = Utility::increment_copy(a, j, -1);

                    auto node_it = nodes_.find(b);
                    if (node_it == std::end(nodes_) or
                        node_it->second.state() != node_state::active)
                        goto add_node_adaptively_skip_iteration;
                }
                forward_neighborhood.insert(a);
                add_node_adaptively_skip_iteration:
            }
        }

        if (forward_neighborhood.empty())
        {
            if (arg & verbosity::empty_forward_neighborhood)
                Print::println("Empty forward neighborhood");

            return error_ret;
        }
        if (arg & verbosity::forward_neighborhood)
        {
            Print::println("Forward_neighborhood:");
            for (const level_type &i: forward_neighborhood)
                Print::println("\t", i);
        }

        fit_nodes(forward_neighborhood);

        if (arg & verbosity::information_of_each_level)
            Print::println("Information of each level:");

        level_type optimal_level {};
        real_t maximum_ratio = 0;
        real_t optimal_level_norm = 0;
        std::unordered_map<level_type, integer_t, Core::hash> coefficients_estimation;
        for (const level_type &level : forward_neighborhood)
        {
            for (const level_type &a : Utility::counter{incremental_})
            {
                if (not std::ranges::equal
                        (level, a, [](index_t arg0, index_t arg1) { return arg0 >= arg1; }))
                    continue;

                level_type al = Utility::subtract_array(level, a);
                index_t parity = std::accumulate(std::cbegin(a), std::cend(a), 0) % 2;
                integer_t c = parity ? -1 : 1;
                coefficients_estimation[al] += c;
            }

            const auto &node = nodes_.find(level)->second;
            real_t norm = node.norm();
            real_t cost = node.cost();
            real_t ratio = norm / cost;
            if (ratio > maximum_ratio)
            {
                maximum_ratio = ratio;
                optimal_level = level;
                optimal_level_norm = norm;
            }

            if (arg & verbosity::information_of_each_level)
                Print::println("\tLevel = ", level, "; Error = ", norm, "; Cost = ", cost,
                               "; ratio = ", ratio);
        }

        if (maximum_ratio == 0)
        {
            Print::println("All nodal contributions are null");
            return {Utility::uniform_array<index_t, n_variates>(negative_1), -1, -1};
        }

        if (arg & verbosity::new_level)
            Print::println("New level = ", optimal_level);

        ASSERT_ASSUME(optimal_level != level_type{});
        activate_node(optimal_level);

        real_t error_estimation = 0;
        for (const auto &i : coefficients_estimation)
            for (const auto &j : coefficients_estimation)
            {
                if (i.first > j.first)
                    continue;
                real_t c = inner_product(i.first, j.first) * i.second * j.second;
                if (i.first == j.first)
                    error_estimation += c;
                else
                    error_estimation += 2 * c;
            }

        ASSERT_ASSUME(Arithmetic::less()(0, error_estimation));
        if (error_estimation < 0)
            error_estimation = 0;
        error_estimation = std::sqrt(error_estimation);

        return {optimal_level, optimal_level_norm, error_estimation};
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    auto smolyak_approximation<function_type, cost_function_type, traits_types...>::
        activate_all_nodes() -> smolyak_approximation &
    {
        std::unordered_set<level_type, Core::hash> levels;

        for (const auto &i : nodes())
            levels.insert(i.first);

        activate_nodes(std::move(levels));

        return *this;
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    template<index_t i, index_t j>
    auto smolyak_approximation<function_type, cost_function_type, traits_types...>::operator_type::
        traits_factory(const level_type &arg0, const traits_initialiser_type &arg1)
        requires(not void_initialiser)
    {
        static_assert(i < n_variates and j < std::tuple_size_v<traits_initialiser_type> and
            (i != 0 or j == 0));

        using traits_type = Utility::get_type<i, traits_types...>;
        index_t level = std::get<i>(arg0);

        if constexpr (std::same_as<typename traits_type::initialiser_type, void>)
        {
            traits_type traits{level};
            if constexpr (i == 0)
                return std::tuple{std::move(traits)};
            else
                return std::tuple_cat(traits_factory<i - 1, j>(arg0, arg1),
                                      std::tuple{std::move(traits)});
        }
        else
        {
            traits_type traits{level, std::get<j>(arg1)};
            if constexpr (i == 0)
                return std::tuple{std::move(traits)};
            else
                return std::tuple_cat(traits_factory<i - 1, j - 1>(arg0, arg1),
                                      std::tuple{std::move(traits)});
        }
    }

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    auto smolyak_approximation<function_type, cost_function_type, traits_types...>::operator_type::
        fit(const function_type &arg0, const cost_function_type &arg1) -> operator_type &
    {
        if (node_state_ != node_state::unfitted)
            return *this;

        auto lambda_coordinates = [](const auto &...traits_secundum)
        {
            return Utility::cartesian_product(traits_secundum.coordinates()...);
        };

        auto lambda_new_coordinates = [](const auto &...traits_secundum)
        {
            return Utility::cartesian_product(traits_secundum.new_coordinates()...);
        };

        std::vector<primary_coordinate_type> coordinates = std::apply(lambda_coordinates, traits_);
        std::vector<primary_coordinate_type> new_coordinates =
            std::apply(lambda_new_coordinates, traits_);

        missing_data missing_data_;

        std::vector<real_t> u, c;
        try
        {
            u = Invocable::evaluate(arg0, coordinates);
        }

        catch (missing_data &e)
        {
            missing_data_.data().merge(e.data());
        }

        try
        {
            c = Invocable::evaluate(arg1, new_coordinates);
        }

        catch (missing_data &e)
        {
            missing_data_.data().merge(e.data());
        }

        if (not std::empty(missing_data_.data()))
            throw missing_data_;

        cost_ = std::accumulate(std::cbegin(c), std::cend(c), 0);
        apply_linear_operator(u, v_, 1);

        node_state_ = node_state::inactive;

        return *this;
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    real_t smolyak_approximation<function_type, cost_function_type, traits_types...>::
        operator_type::evaluate() const requires(void_secondary_domain)
    {
        ASSERT_ASSUME(node_state_ != node_state::unfitted);
        std::vector<real_t> v = v_, w;
        evaluate(v, w);
        ASSERT_ASSUME(std::size(w) == 1);

        return w.front();
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    real_t smolyak_approximation<function_type, cost_function_type, traits_types...>::
        operator_type::evaluate(const secondary_coordinate_type &x) const
        requires(not void_secondary_domain)
    {
        ASSERT_ASSUME(node_state_ != node_state::unfitted);
        std::vector<real_t> v = v_, w;
        evaluate(v, w, x);
        ASSERT_ASSUME(std::size(w) == 1);

        return w.front();
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    std::vector<real_t> smolyak_approximation<function_type, cost_function_type, traits_types...>::
        operator_type::homomorphism(const homomorphism_c auto &...homomorphisms) const
    requires(sizeof...(homomorphisms) == n_variates)
    {
        ASSERT_ASSUME(node_state_ != node_state::unfitted);
        std::vector<real_t> v = v_, w;
        index_t max_size = homomorphism_apparatus_maximum_size(std::size(v), 0, homomorphisms...);
        v.reserve(max_size);
        w.reserve(max_size);
        homomorphism(v, w, 1, homomorphisms...);

        return w;
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
        template<index_t i_variate>
    void smolyak_approximation<function_type, cost_function_type, traits_types...>::operator_type::
    apply_linear_operator(std::vector<real_t> &u, std::vector<real_t> &v, index_t stride)
    {
        using current_traits_type = Utility::get_type<i_variate, traits_types...>;
        const current_traits_type &current_traits = std::get<i_variate>(traits_);
        index_t s = std::size(current_traits.coordinates());

        if constexpr (Apparatus::smolyak_traits_trivial_linear_operator_c<current_traits_type>)
        {
            if constexpr (i_variate != 0)
                apply_linear_operator<i_variate - 1>(u, v, stride * s);

            else
                std::swap(u, v);

            return;
        }

        else
        {
            index_t t = current_traits.linear_operator_output_size();

            ASSERT_ASSUME(std::size(u) % s == 0);
            index_t output_size = std::size(u) / s * t;

            v.resize(output_size);

            for (auto it_u = std::begin(u), it_v = std::begin(v);
                 it_u < std::end(u); it_u += s * stride, it_v += t * stride)
                for (auto it_u_secundum = it_u, it_v_secundum = it_v;
                     it_u_secundum != it_u + stride;
                     ++it_u_secundum, ++it_v_secundum)
                    current_traits.linear_operator
                        (stride,
                            std::span<const real_t>{it_u_secundum,
                                                    it_u_secundum + (s - 1) * stride + 1},
                            std::span<real_t>{it_v_secundum, it_v_secundum + (t - 1) * stride + 1});

            if constexpr (i_variate != 0)
            {
                std::swap(u, v);
                apply_linear_operator<i_variate - 1>(u, v, stride * t);
            }

            return;
        }
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
        template<index_t i_variate>
    void smolyak_approximation<function_type, cost_function_type, traits_types...>::operator_type::
    evaluate(std::vector<real_t> &v, std::vector<real_t> &w) const
    requires(void_secondary_domain)
    {
        using current_traits_type = Utility::get_type<i_variate, traits_types...>;
        const current_traits_type &current_traits = std::get<i_variate>(traits_);

        index_t t;
        if constexpr (Apparatus::smolyak_traits_trivial_linear_operator_c<current_traits_type>)
            t = std::size(current_traits.coordinates());
        else
            t = current_traits.linear_operator_output_size();

        ASSERT_ASSUME(std::size(v) % t == 0);
        w.resize(std::size(v) / t);

        for (auto it_v = std::begin(v), it_w = std::begin(w); it_w != std::end(w);
             it_v += t, ++it_w)
            *it_w = current_traits.evaluate(std::span<const real_t>{it_v, it_v + t});

        if constexpr (i_variate != 0)
        {
            std::swap(v, w);
            evaluate<i_variate - 1>(v, w);
        }

        return;
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
        template<index_t i_variate, index_t j_variate>
    void smolyak_approximation<function_type, cost_function_type, traits_types...>::operator_type::
        evaluate(std::vector<real_t> &v, std::vector<real_t> &w,
        const secondary_coordinate_type &x) const
    requires(not void_secondary_domain)
    {
        using current_traits_type = Utility::get_type<i_variate, traits_types...>;
        const current_traits_type &current_traits = std::get<i_variate>(traits_);

        constexpr bool current_void_secondary_domain =
                std::same_as<typename current_traits_type::secondary_domain, void>;

        index_t t;
        if constexpr (Apparatus::smolyak_traits_trivial_linear_operator_c<current_traits_type>)
            t = std::size(current_traits.coordinates());
        else
            t = current_traits.linear_operator_output_size();
        ASSERT_ASSUME(std::size(v) % t == 0);
        w.resize(std::size(v) / t);

        if constexpr (current_void_secondary_domain)
            for (auto it_v = std::begin(v), it_w = std::begin(w); it_w != std::end(w);
                 it_v += t, ++it_w)
                *it_w = current_traits.evaluate(std::span<const real_t>{it_v, it_v + t});

        else
        {
            using current_secondary_domain_variate_type =
                    std::tuple_element_t<j_variate, secondary_coordinate_type>;
            const current_secondary_domain_variate_type &current_secondary_domain_variate =
                    std::get<j_variate>(x);

            for (auto it_v = std::begin(v), it_w = std::begin(w); it_w != std::end(w);
                    it_v += t, ++it_w)
                *it_w = current_traits.evaluate
                (std::span<const real_t>{it_v, it_v + t}, current_secondary_domain_variate);
        }

        if constexpr (i_variate != 0)
        {
            std::swap(v, w);
            evaluate<i_variate - 1, j_variate - not current_void_secondary_domain>(v, w, x);
        }

        return;
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    template<index_t i_variate>
    index_t smolyak_approximation<function_type, cost_function_type, traits_types...>::
        operator_type:: homomorphism_apparatus_maximum_size(index_t size, index_t max_size,
            const homomorphism_c auto &...homomorphisms) const
            requires(sizeof...(homomorphisms) == n_variates)
    {
        const auto &homomorphism_ = Utility::get_value<i_variate>(homomorphisms...);
        index_t level = std::get<i_variate>(level_);

        index_t s = homomorphism_.input_size(level);
        index_t t = homomorphism_.output_size(level);
        index_t new_size = size / s * t;
        max_size = std::max(new_size, max_size);

        if constexpr (i_variate != 0)
            max_size = homomorphism_apparatus_maximum_size<i_variate - 1>
                (new_size, max_size, homomorphisms...);

        return max_size;
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    template<index_t i_variate>
    void smolyak_approximation<function_type, cost_function_type, traits_types...>::operator_type::
        homomorphism(std::vector<real_t> &v, std::vector<real_t> &w, index_t stride,
            const homomorphism_c auto &...homomorphisms) const
            requires(sizeof...(homomorphisms) == n_variates)
    {
        const auto &homomorphism_ = Utility::get_value<i_variate>(homomorphisms...);
        index_t level = std::get<i_variate>(level_);

        index_t s = homomorphism_.input_size(level);
        index_t t = homomorphism_.output_size(level);
        ASSERT_ASSUME(std::size(v) % s == 0);
        ASSERT_ASSUME(w.capacity() >= std::size(v) / s * t);
        w.resize(std::size(v) / s * t);

        for (auto it_v = std::begin(v), it_w = std::begin(w);
             it_w < std::end(w);
             it_v += s * stride, it_w += t * stride)
            for (auto it_v_secundum = it_v, it_w_secundum = it_w;
                 it_v_secundum != it_v + stride;
                 ++it_v_secundum, ++it_w_secundum)
                homomorphism_
                        (level, stride,
                         std::span<const real_t>{it_v_secundum,
                                                 it_v_secundum + (s - 1) * stride + 1},
                         std::span<real_t>{it_w_secundum, it_w_secundum + (t - 1) * stride + 1});

        if constexpr (i_variate != 0)
        {
            std::swap(v, w);
            stride *= t;
            homomorphism<i_variate - 1>(v, w, stride, homomorphisms...);
        }

        return;
    }

    template<class function_type, class cost_function_type, smolyak_traits_c... traits_types>
    real_t smolyak_approximation<function_type, cost_function_type, traits_types...>::
        operator_type::inner_product(const operator_type &arg) const
    {
        static auto lambda_embedding = [](const auto &traits_secundum, index_t level)
        {
            return traits_secundum.embedding(level);
        };

        static auto lambda_dual_embedding = [](const auto &traits_secundum, index_t level)
        {
            return traits_secundum.dual_embedding(level);
        };

        static auto lambda_homomorphism = [](const operator_type &node)
        {
            return [&node](const auto &...args)
            {
                return node.homomorphism(args...);
            };
        };

        level_type common_level = Utility::max_array(level_, arg.level_);

        std::tuple embeddings{Utility::transform_tuple(lambda_embedding, traits_, common_level)};
        std::vector<real_t> vector = std::apply(lambda_homomorphism(*this), embeddings);

        std::tuple dual_embeddings{Utility::transform_tuple
            (lambda_dual_embedding, arg.traits_, common_level)};
        std::vector<real_t> dual_vector = std::apply(lambda_homomorphism(arg), dual_embeddings);

        ASSERT_ASSUME(std::size(vector) == std::size(dual_vector));

        return std::inner_product
            (std::cbegin(vector), std::cend(vector), std::cbegin(dual_vector),
             static_cast<real_t>(0));
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

