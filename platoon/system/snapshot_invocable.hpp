#pragma once

#include "fwd.hpp"
#include "invocable/invocable_disk_binary_tree.hpp"

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Platoon
{
    std::string file_name(const std::array<index_t, 5> &arg)
    {
        return std::string{"../data/"} + std::to_string(arg[0]) + "_" + std::to_string(arg[1]) +
               "_" + std::to_string(arg[2]) + "_" + std::to_string(arg[3]) + "_" +
               std::to_string(arg[4]);
    }

    std::string file_name_data(const std::array<index_t, 5> &arg)
    {
        return file_name(arg) + "/data";
    }

    bool exists(const std::array<index_t, 5> &arg)
    {
        return std::filesystem::exists(file_name(arg));
    }

    bool exists_data(const std::array<index_t, 5> &arg)
    {
        return std::filesystem::exists(file_name_data(arg));
    }

    bool exists_data(const input_type &arg)
    {
        std::array<index_t, 5> arg_bis;
        arg_bis[0] = std::get<1>(arg);
        arg_bis[1] = std::get<2>(arg);
        arg_bis[2] = sampling_type::index(std::get<3>(arg));
        arg_bis[3] = sampling_type::index(std::get<4>(arg));
        arg_bis[4] = sampling_type::index(std::get<5>(arg));

        return exists_data(arg_bis);
    }

    class snapshot_data
    {
    private:
        index_t component_;
        index_t fidelity_;
        index_t geometry_;
        index_t spacing_index_;
        index_t yaw_index_;
        index_t reynolds_index_;
        real_t spacing_;
        real_t yaw_;
        real_t reynolds_;
        real_t qoi_;
        real_t cost_;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        snapshot_data() = delete;

        snapshot_data(index_t component, index_t fidelity, index_t geometry, index_t spacing_index,
            index_t yaw_index, index_t reynolds_index, real_t qoi, real_t cost) :
            component_{component},
            fidelity_{fidelity},
            geometry_{geometry},
            spacing_index_{spacing_index},
            yaw_index_{yaw_index},
            reynolds_index_{reynolds_index},
            spacing_{sampling_type::coordinate_from_index(spacing_index_)},
            yaw_{sampling_type::coordinate_from_index(yaw_index_)},
            reynolds_{sampling_type::coordinate_from_index(reynolds_index_)},
            qoi_{qoi},
            cost_{cost} {}

        operator input_type() const
        {
            return{component_, fidelity_, geometry_, spacing_, yaw_, reynolds_};
        }

        operator output_type() const
        {
            return std::array<real_t, 2>{qoi_, cost_};
        }

        operator pair_type() const
        {
            return {static_cast<input_type>(*this), static_cast<output_type>(*this)};
        }

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        index_t &component() { return component_; }
        const index_t &component() const { return component_; }
        index_t &fidelity() { return fidelity_; }
        const index_t &fidelity() const { return fidelity_; }
        index_t &geometry() { return geometry_; }
        const index_t &geometry() const { return geometry_; }
        index_t &spacing_index() { return spacing_index_; }
        const index_t &spacing_index() const { return spacing_index_; }
        index_t &yaw_index() { return yaw_index_; }
        const index_t &yaw_index() const { return yaw_index_; }
        index_t &reynolds_index() { return reynolds_index_; }
        const index_t &reynolds_index() const { return reynolds_index_; }
        real_t &spacing() { return spacing_; }
        const real_t &spacing() const { return spacing_; }
        real_t &yaw() { return yaw_; }
        const real_t &yaw() const { return yaw_; }
        real_t &reynolds() { return reynolds_; }
        const real_t &reynolds() const { return reynolds_; }
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
    void read<Platoon::snapshot_data>
        (std::ifstream &file, Platoon::snapshot_data &value)
    {
        Arithmetic::index_t component, fidelity, geometry, spacing_index, yaw_index, reynolds_index;
        double qoi, cost;

        Disk::read(file, component);
        Disk::read(file, fidelity);
        Disk::read(file, geometry);
        Disk::read(file, spacing_index);
        Disk::read(file, yaw_index);
        Disk::read(file, reynolds_index);
        Disk::read(file, qoi);
        Disk::read(file, cost);

        value = Platoon::snapshot_data
            {component, fidelity, geometry, spacing_index, yaw_index, reynolds_index,
             static_cast<Arithmetic::real_t>(qoi), static_cast<Arithmetic::real_t>(cost)};

        return;
    }

    template<>
    void write<Platoon::snapshot_data>
        (std::ofstream &file, const Platoon::snapshot_data &value)
    {
        Disk::write(file, value.component());
        Disk::write(file, value.fidelity());
        Disk::write(file, value.geometry());
        Disk::write(file, value.spacing_index());
        Disk::write(file, value.yaw_index());
        Disk::write(file, value.reynolds_index());
        Disk::write(file, static_cast<double>(value.qoi()));
        Disk::write(file, static_cast<double>(value.cost()));

        return;
    }
}

// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
// **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //

namespace Platoon
{
    class snapshot_invocable
    {
    private:
        static constexpr std::array<std::array<real_t, 2>, 3> normalisation_factors_
        {std::array<real_t, 2>{0.142089, 0.155349}, std::array<real_t, 2>{0.210775, -0.00213132}, std::array<real_t, 2>{0.110602, -0.0666261}};
    public:
        using invocable_type = Invocable::invocable_disk_binary_tree
            <std::array<double, 4>(index_t, index_t, double, double, double)>;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    private:
        invocable_type binary_tree_{"data"};

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        snapshot_invocable() = default;
        ~snapshot_invocable() = default;

        snapshot_invocable(const snapshot_invocable &) = delete;
        snapshot_invocable &operator=(const snapshot_invocable &) = delete;

        snapshot_invocable(snapshot_invocable &&) = delete;
        snapshot_invocable &operator=(snapshot_invocable &&) = delete;

        // **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** **** //
    public:
        static const auto &normalisation_factors()
        {
            return normalisation_factors_;
        }

        std::array<real_t, 4> operator()(index_t fidelity, index_t geometry, real_t spacing,
            real_t yaw, real_t reynolds)
        {
            std::array<double, 4> array = binary_tree_(fidelity, geometry,
                static_cast<double>(spacing), static_cast<double>(yaw),
                static_cast<double>(reynolds));

            std::array<real_t, 4> ret;

            std::copy(std::cbegin(array), std::cend(array), std::begin(ret));

            return ret;
        }

        auto qoi_function()
        {
            return [this](const input_type &arg) -> real_t
            {
                index_t c = std::get<0>(arg);

                return ((*this)(std::get<1>(arg), std::get<2>(arg), std::get<3>(arg),
                    std::get<4>(arg), std::get<5>(arg))[c] -
                        normalisation_factors()[c][1]) / normalisation_factors()[c][0];
            };
        }

        auto cost_function()
        {
            return [this](const input_type &arg) -> real_t
            {
                return (*this)(std::get<1>(arg), std::get<2>(arg), std::get<3>(arg),
                               std::get<4>(arg), std::get<5>(arg)).back();
            };
        }
    };
}

