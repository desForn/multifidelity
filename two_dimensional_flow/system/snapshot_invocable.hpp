#pragma once

#include "fwd.hpp"
#include "invocable/invocable_disk_binary_tree.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Two_dimensional_flow
{
    template<class sampling_type>
    class snapshot_data
    {
    private:
        index_t fidelity_;
        real_t angle_;
        index_t angle_index_;
        real_t qoi_;
        real_t cost_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        snapshot_data() = default;

        snapshot_data(index_t fidelity, real_t angle) :
                fidelity_{fidelity},
                angle_{angle},
                angle_index_{sampling_type::index(angle_)} {}

        snapshot_data(index_t fidelity, real_t angle, real_t qoi, real_t cost) :
                fidelity_{fidelity},
                angle_{angle},
                angle_index_{sampling_type::index(angle_)},
                qoi_{qoi},
                cost_{cost} {}

        operator fidelity_coordinate_type() const
        {
            return {std::array{fidelity_}, std::array{angle_}};
        }

        operator fidelity_coordinate_index_type() const
        {
            return {std::array{fidelity_}, std::array{angle_index_}};
        }

        operator pair_type() const
        {
            return {static_cast<fidelity_coordinate_type>(*this), output_type{qoi_, cost_}};
        }

        operator std::string() const
        {
            return std::to_string(fidelity_) + "_" + std::to_string(angle_index_);
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        index_t &fidelity() { return fidelity_; }
        const index_t &fidelity() const { return fidelity_; }

        real_t &angle() { return angle_; }
        const real_t &angle() const { return angle_; }

        index_t &angle_index() { return angle_index_; }
        const index_t &angle_index() const { return angle_index_; }

        real_t &qoi() { return qoi_; }
        const real_t &qoi() const { return qoi_; }

        real_t &cost() { return cost_; }
        const real_t &cost() const { return cost_; }
    };
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Disk
{
    template<>
    void read<Two_dimensional_flow::pair_type>(std::ifstream &file, Two_dimensional_flow::pair_type &value)
    {
        Disk::read(file, std::get<0>(std::get<0>(std::get<0>(value))));
        Disk::read(file, std::get<0>(std::get<1>(std::get<0>(value))));
        Disk::read(file, std::get<0>(std::get<1>(value)));
        Disk::read(file, std::get<1>(std::get<1>(value)));

        return;
    }

    template<>
    void write<Two_dimensional_flow::pair_type>(std::ofstream &file, const Two_dimensional_flow::pair_type &value)
    {
        Disk::write(file, std::get<0>(std::get<0>(std::get<0>(value))));
        Disk::write(file, std::get<0>(std::get<1>(std::get<0>(value))));
        Disk::write(file, std::get<0>(std::get<1>(value)));
        Disk::write(file, std::get<1>(std::get<1>(value)));

        return;
    }

    template<>
    void read<Two_dimensional_flow::snapshot_data<Sampling::exponential<Sampling::chebyshev_lobatto>>>
            (std::ifstream &file,
             Two_dimensional_flow::snapshot_data<Sampling::exponential<Sampling::chebyshev_lobatto>> &value)
    {
        Disk::read(file, value.fidelity());
        Disk::read(file, value.angle());
        value.angle_index() =
                Sampling::exponential<Sampling::chebyshev_lobatto>::index(value.angle());
        Disk::read(file, value.qoi());
        Disk::read(file, value.cost());

        return;
    }

    template<>
    void write<Two_dimensional_flow::snapshot_data<Sampling::exponential<Sampling::chebyshev_lobatto>>>
            (std::ofstream &file,
             const Two_dimensional_flow::snapshot_data<Sampling::exponential<Sampling::chebyshev_lobatto>> &
             value)
    {
        Disk::write(file, value.fidelity());
        Disk::write(file, value.angle());
        Disk::write(file, value.qoi());
        Disk::write(file, value.cost());

        return;
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Two_dimensional_flow
{
    class snapshot_invocable
    {
    public:
        using fidelity_type = std::array<index_t, 1>;
        using coordinate_type = std::array<real_t, 1>;
        using fidelity_coordinate_type = std::tuple<fidelity_type, coordinate_type>;

        using output_type = std::array<real_t, 2>;

        using invocable_type = Invocable::invocable_disk_binary_tree
            <output_type(fidelity_type, coordinate_type)>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        invocable_type binary_tree_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        snapshot_invocable() = delete;
        ~snapshot_invocable() = default;

        snapshot_invocable(const snapshot_invocable &) = delete;
        snapshot_invocable &operator=(const snapshot_invocable &) = delete;

        snapshot_invocable(snapshot_invocable &&) = delete;
        snapshot_invocable &operator=(snapshot_invocable &&) = delete;

        snapshot_invocable(const std::string &);

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        output_type operator()(const fidelity_type &, const coordinate_type &);

        auto qoi_function();
        auto cost_function();
    };

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    snapshot_invocable::snapshot_invocable(const std::string &data_file) :
            binary_tree_{data_file} {}

    // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

    auto snapshot_invocable::operator()(const fidelity_type &arg0, const coordinate_type &arg1)
    -> output_type
    {
        return binary_tree_(arg0, arg1);
    }

    auto snapshot_invocable::qoi_function()
    {
        return [this](const fidelity_type &arg0, const coordinate_type &arg1)
        {
            return std::get<0>((*this)(arg0, arg1));
        };
    }

    auto snapshot_invocable::cost_function()
    {
        return [this](const fidelity_type &arg0, const coordinate_type &arg1)
        {
            return std::get<1>((*this)(arg0, arg1));
        };
    }
}
