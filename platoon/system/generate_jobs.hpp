#pragma once

#include "snapshot_invocable.hpp"

namespace Platoon
{
    void generate_jobs(const std::vector<input_type> &arg)
    {
        std::string str;
        const char *c_str;

        auto sort = []
            (const std::array<index_t, Domain::n_parameters> &arg0,
             const std::array<index_t, Domain::n_parameters> &arg1) -> bool
            {
                index_t sum0 = std::accumulate(std::cbegin(arg0), std::cend(arg0), 0);
                index_t sum1 = std::accumulate(std::cbegin(arg1), std::cend(arg1), 0);

                if (sum0 < sum1)
                    return true;

                if (sum0 > sum1)
                    return false;

                else
                    return arg0 < arg1;
            };

        std::array<std::set<std::array<index_t, Domain::n_parameters>, decltype(sort)>,
            Domain::n_fidelity_levels> sets;

        for (const input_type &i : arg)
            sets[std::get<1>(i)].insert(std::array<index_t, 4>
                {std::get<2>(i), sampling_type::index(std::get<3>(i)),
                 sampling_type::index(std::get<4>(i)),
                 sampling_type::index(std::get<5>(i))});

        std::array<std::vector<std::ofstream>, Domain::n_fidelity_levels> files;

        for (index_t i = 0; i != Domain::n_fidelity_levels; ++i)
            for (index_t j = 0; j != std::min(Settings::max_jobs[i], std::size(sets[i])); ++j)
            {
                ASSERT_ASSUME(Settings::n_processors[i] != 0);

                files[i].emplace_back("job" + std::to_string(i) + "_" + std::to_string(j));
                ASSERT_ASSUME(files[i].back().is_open());

                str = "#!/bin/bash\n";
                str += "#SBATCH --job-name=" + std::to_string(i) + "_" + std::to_string(j) + "\n";

                if (Settings::partition[i] == 1)
                {
                    str += "#SBATCH --partition=standard\n";
                    str += "#SBATCH --ntasks=" + std::to_string(Settings::n_processors[i]) + "\n";
                }

                else if (Settings::partition[i] == 16 or Settings::partition[i] == 28 or
                    Settings::partition[i] == 64)
                {
                    ASSERT_ASSUME(Settings::n_processors[i] % Settings::partition[i] == 0);
                    str += "#SBATCH --partition=mpi" + std::to_string(Settings::partition[i]) +
                           "c\n";
                    str += "#SBATCH --ntasks=" + std::to_string(Settings::n_processors[i]) +
                           "\n";
                    str += "#SBATCH --nodes=" +
                        std::to_string(Settings::n_processors[i] / Settings::partition[i]) + "\n";
                }

                else
                    ASSERT_ASSUME(false);

                str += "echo \"Start:\" $(date)\n\n";

                c_str = str.c_str();

                files[i].back().write(c_str, std::strlen(c_str));
            }

        for (index_t i = 0; i != Domain::n_fidelity_levels; ++i)
        {
            index_t j = 0;
            bool parallel = Settings::n_processors[i] != 1;

            for (const std::array<index_t, Domain::n_parameters> &k : sets[i])
            {
                real_t spacing = sampling_type::coordinate_from_index(k[1]);
                spacing = (spacing + 1) / 2;
                spacing =
                    Domain::min_spacing + (Domain::max_spacing - Domain::min_spacing) * spacing;

                real_t yaw = sampling_type::coordinate_from_index(k[2]);
                yaw = (yaw + 1) / 2;
                yaw = Domain::min_yaw + (Domain::max_yaw - Domain::min_yaw) * yaw;

                real_t reynolds = sampling_type::coordinate_from_index(k[3]);
                reynolds = (reynolds + 1) / 2;
                reynolds =
                    Domain::min_reynolds + (Domain::max_reynolds - Domain::min_reynolds) * reynolds;

                std::string file_name = std::string{"../data/"} + std::to_string(i);
                for (index_t a: k)
                    file_name += "_" + std::to_string(a);

                if (std::filesystem::exists(file_name))
                {
                    Print::println(file_name, "\nalready exists.");
                    std::abort();
                }

                str = "begin=$(date +%s)\n";
                str += "cp -r meshes/s" + std::to_string(k.front()) + "f" +
                    std::to_string(i) + " " + file_name + "\n";
                str += "cd " + file_name + "\n";

                str += "../../system/generate_bcs " + std::to_string(yaw) + "\n";
                str += "../../system/generate_nu " + std::to_string(reynolds) + "\n";

                if (parallel)
                {
                    str += "mpiexec ../../system/mesh_displacement/mesh_displacement -s " +
                        std::to_string(spacing) + " -g " + std::to_string(k.front()) +
                        " -parallel > log.mesh_displacement\n";
                    str += "for i in {0.." + std::to_string(Settings::n_processors[i] - 1) + "}\n";
                    str += "do\n";
                    str += "    rm processor${i}/constant/polyMesh/points\n";
                    str += "    mv processor${i}/0/polyMesh/points";
                    str += " processor${i}/constant/polyMesh/points\n";
                    str += "    rm -r processor${i}/0/polyMesh\n";
                    str += "done\n";
                    if (i != 0)
                    {
                        std::array<index_t, Domain::n_parameters + 1> previous_fidelity;
                        previous_fidelity.front() = i - 1;

                        std::copy(std::cbegin(k), std::cend(k), std::begin(previous_fidelity) + 1);

                        str += "mapFields " + Platoon::file_name(previous_fidelity);
                        str += " -consistent -sourceTime latestTime -parallelTarget";
                        if (Settings::n_processors[i - 1] != 1)
                            str += " -parallelSource";
                        str += " > log.mapFields\n";
                    }
                    str += "mpiexec simpleFoam -parallel > log.simpleFoam\n";
                }

                else
                {
                    str += "../../system/mesh_displacement/mesh_displacement -s " +
                           std::to_string(spacing) + " -g " + std::to_string(k.front()) +
                           " > log.mesh_displacement\n";
                    str += "rm constant/polyMesh/points\n";
                    str += "mv 0/polyMesh/points constant/polyMesh/points\n";
                    str += "rm -r 0/polyMesh\n";
                    str += "simpleFoam > log.simpleFoam\n";
                }

                str += "end=$(date +%s)\n";
                str += "data=$(../../system/qoi " + std::to_string(i);
                str += " \"$(((end - begin) * " + std::to_string(Settings::n_processors[i]);
                str += "))\")\n";
                str += "echo \"$data\" >> data_ascii\n\n";

                c_str = str.c_str();
                files[i][j].write(c_str, std::strlen(c_str));

                ++j;

                if (j == Settings::max_jobs[i])
                    j = 0;
            }
        }

        for (auto &i : files)
            for (auto &j : i)
            {
                str = "echo \"End    : $(date)\"\n\n";
                c_str = str.c_str();
                j.write(c_str, std::strlen(c_str));

                j.close();
            }

        return;
    }
}

