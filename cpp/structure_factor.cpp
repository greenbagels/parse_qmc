/*! @file structure_factor.cpp
 *  @brief Calculates structure factors from averaged QMC correlators
 *  @author Sameed Pervaiz (pervaiz.8@osu.edu)
 *  @date 2021-09-19
 *  @copyright GPLv3
 */

#include <iostream>
#include <filesystem>
#include <fstream>
#include <map>
#include <cmath>
#include <sstream>
#include "common.hpp"
#include "quantity.hpp"

using namespace qmc;

enum fixed_qtys
{
    MU = 0,
    BETA,
    NUM_FIXED
};

enum measured_qtys
{
    DENSITY = 0,
    DOUBLON,
    TOTAL_SIGN,
    NUM_MEASURED
};

enum cf_fields
{
    MOMENT = 0,
    SPIN,
    NUM_CFS
};

// this represents data from a single run.
struct qmc_data_full
{
    qmc_data_full(size_t n) : cfs{matrix<double>(n), matrix<double>(n)} {}
    quantity<double> fixed[NUM_FIXED];
    quantity<double> measured[NUM_MEASURED];
    matrix<double> cfs[NUM_CFS];
};

// T requirements: has .value(), .quantity() of numeric type
// template <>
qmc_data_full average(std::vector<qmc_data_full> &list)
{
    qmc_data_full data(list[0].cfs[0].width());
    // First, the data that is zero-uncertainty and constant across all entries
    data.fixed[MU] = list[0].fixed[MU];
    data.fixed[BETA] = list[0].fixed[BETA];
    // Now, do the averaging over the measured entries
    // The ordering will thrash cache, but lets us automatically stop when we
    // encounter a zero uncertainty value.
    double weighted_sum(0.), sum_of_weights(0.);
    for (int i = 0; i < NUM_MEASURED; i++)
    {
        double weighted_sum(0.), sum_of_weights(0.);
        for (auto &qdata : list)
        {
            // Take zero-uncertainty at face value and move on to the next field
            if (qdata.measured[i].uncertainty() == 0.)
            {
                data.measured[i] = qdata.measured[i];
                goto CONTINUE_MEASURED;
            }
            auto weight = 1. / (qdata.measured[i].uncertainty() * qdata.measured[i].uncertainty());
            weighted_sum += weight * qdata.measured[i].value();
            sum_of_weights += weight;
        }
        data.measured[i] = quantity<double>(weighted_sum / sum_of_weights, 1. / std::sqrt(sum_of_weights));
        CONTINUE_MEASURED:
        continue;
    }
    // Now, the CFs
    for (int i = 0; i < NUM_CFS; i++)
    {
        for (auto y = 0; y < list[0].cfs[i].width(); y++)
        {
            for (auto x = 0; x < list[0].cfs[i].width(); x++)
            {
                double weighted_sum(0.), sum_of_weights(0.);
                for (auto &qdata : list)
                {
                    if (qdata.cfs[i](x,y).uncertainty() == 0.)
                    {
                        data.cfs[i](x,y) = qdata.cfs[i](x,y);
                        goto CONTINUE_CFS;
                    }
                    auto weight = 1. / (qdata.cfs[i](x,y).uncertainty() * qdata.cfs[i](x,y).uncertainty());
                    weighted_sum += weight * qdata.cfs[i](x,y).value();
                    sum_of_weights += weight;
                }
                data.cfs[i](x,y) = quantity<double>(weighted_sum / sum_of_weights, 1. / std::sqrt(sum_of_weights));
                CONTINUE_CFS:
                continue;
            }
        }
    }
    return data;
}

template <typename T>
quantity<T> parse_line_with_uc(std::ifstream &file, std::string field)
{
    for (std::string line; std::getline(file, line);)
    {
        if (line.find(field) != std::string::npos)
        {
            T val, uc;
            if (line.find("=") == std::string::npos)
                throw std::runtime_error("Ill-formed data line in file:\n" + line);
            line = line.substr(line.find("=")+1);
            std::stringstream sstr(line);
            sstr >> val >> uc;
            return quantity<T>(val, uc);
        }
    }
    throw std::runtime_error("Error: could not find line-data field name in parsed data file!");
}

// TODO: "smartly" parameterize n
template <typename T>
matrix<T> parse_matrix_with_uc(std::ifstream &file, std::string field, size_t n)
{
    for (std::string line; std::getline(file, line);)
    {
        if (line.find(field) != std::string::npos)
        {
            matrix<T> mat(n);
            for (auto i = n*n; i > 0; i--)
            {
                std::getline(file, line);
                std::stringstream sstr(line);
                size_t x,y;
                T val, uc;
                sstr >> x >> y >> val >> uc;
                mat(x,y) = quantity<T>(val, uc);
            }
            return mat;
        }
    }
    throw std::runtime_error("Error: could not find matrix-data field name in parsed data file!");
}

template <typename T>
int parse_qmc_file_full(T file_handle, qmc_data_full &data)
{
    std::ifstream qmc_file(file_handle);
    if (!qmc_file.is_open())
    {
        std::cerr << "Error: could not open file\n";
        return 1;
    }
    data.measured[TOTAL_SIGN] = parse_line_with_uc<double>(qmc_file, "Average total sign");
    data.measured[DENSITY] = parse_line_with_uc<double>(qmc_file, "Average density");
    data.measured[DOUBLON] = parse_line_with_uc<double>(qmc_file, "Average Nup*Ndn");
    data.cfs[MOMENT] = parse_matrix_with_uc<double>(qmc_file, "mi2x-mi2x correlation function", data.cfs[MOMENT].width());
    data.cfs[SPIN] = parse_matrix_with_uc<double>(qmc_file, "zz Spin correlation function", data.cfs[SPIN].width());
    return 0;
}

matrix<double> calc_structure_factors(qmc_data_full &data, int type, int n)
{
    auto width = data.cfs[type].width();
    matrix<double> sfdata(n);
    const auto pi = std::acos(-1);
    for (auto iqy = 0; iqy < n; iqy++)
    {
        for (auto iqx = 0; iqx < n; iqx++)
        {
            double qx = iqx * 2. * pi / n;
            double qy = iqy * 2. * pi / n;
            auto total = 0;
            for (auto iy = 0; iy < width; iy++)
            {
                for (auto ix = 0; ix < width; ix++)
                {
                    for (auto jy = 0; jy < width; jy++)
                    {
                        for (auto jx = 0; jx < width; jx++)
                        {
                            auto rx = ix-jx, ry = iy-jy;
                            auto absrx = std::abs(ix-jx), absry = std::abs(iy-jy);
                            auto num_pairs = (n - absrx) * (n - absry);
                            total += num_pairs;
                            auto cf = data.cfs[type](absrx, absry);
                            sfdata(iqx, iqy) += num_pairs * std::cos(rx * qx + ry * qy) * cf;
                        }
                    }
                }
            }
            sfdata(iqx, iqy) /= total;
        }
    }
    return sfdata;
}

int main(int argc, char* argv[])
{
    namespace fs = std::filesystem;
    // This is the loop over beta
    std::ofstream max_spin_sf_file("max_spin_sf.dat");
    std::ofstream max_moment_sf_file("max_moment_sf.dat");
    for (auto &p : fs::directory_iterator(argv[1]))
    {
        double beta;
        quantity<double> max_spin_sf, max_moment_sf;
        if (!p.is_directory())
        {
            std::cerr << "Skipping non-directory file " << p.path().filename() << std::endl;
            continue;
        }
        fs::current_path(p);
        // This is the loop over mu and r.
        std::map<double, std::vector<qmc_data_full>> unaveraged_data;
        for (auto &p : fs::directory_iterator(fs::current_path()))
        {
            auto fname = p.path().filename().string();
            if (fname.find("rz") == std::string::npos)
            {
                //std::cerr << "Warning: skipping non-data file " << fname << std::endl;
                continue;
            }
            qmc_data header;
            try
            {
                header = parse_qmc_fname(fname);
            }
            catch (std::exception &e)
            {
                std::cerr << "Exception caught in filename parsing: " <<  e.what() << std::endl;
                return 1;
            }
            unaveraged_data[header.mu].emplace_back(header.n/2 + 1);
            auto &data = unaveraged_data[header.mu].back();
            data.fixed[MU] = header.mu;
            data.fixed[BETA] = header.beta;
            beta = header.beta;
            parse_qmc_file_full(p.path(), data);
        }
        // The files are parsed, so now average and calc sfs.
        std::map<double, qmc_data_full> averaged_data;
        std::ofstream moment_sf_file("../../moment_sf_beta_" + std::to_string(beta) + ".dat");
        std::ofstream spin_sf_file("../../spin_sf_beta_" + std::to_string(beta) + ".dat");
        for (auto i = unaveraged_data.begin(); i != unaveraged_data.end(); i++)
        {
            auto mu = i->first;
            auto &data = i->second;
            averaged_data.emplace(mu, ::average(data));
            // std::cerr << averaged_data.at(mu).measured[DENSITY] << std::endl;
            auto moment_sf = calc_structure_factors(averaged_data.at(mu), MOMENT, 10);
            auto spin_sf = calc_structure_factors(averaged_data.at(mu), SPIN, 10);
            for (auto iqy = 0; iqy < moment_sf.width(); iqy++)
            {
                for (auto iqx = 0; iqx < moment_sf.width(); iqx++)
                {
                    const auto qy = iqy * 2 * std::acos(-1) / moment_sf.width();
                    const auto qx = iqx * 2 * std::acos(-1) / moment_sf.width();
                    moment_sf_file << qy << " " << qx << " " << mu << " "
                    << moment_sf(iqx,iqy) << "\n";
                    spin_sf_file << qy << " " << qx << " " << mu << " "
                    << spin_sf(iqx,iqy) << "\n";
                    max_spin_sf = std::max(spin_sf(iqx, iqy), max_spin_sf,
                            [](quantity<double> a, quantity<double> b){return a.value() < b.value();});
                    max_moment_sf = std::max(moment_sf(iqx, iqy), max_moment_sf,
                            [](quantity<double> a, quantity<double> b){return a.value() < b.value();});
                }
            }
            moment_sf_file << "\n";
            spin_sf_file << "\n";
        }
        max_spin_sf_file << beta << " " << max_spin_sf << std::endl;
        max_moment_sf_file << beta << " " << max_moment_sf << std::endl;
        fs::current_path(fs::current_path().parent_path().parent_path());
    }
    return 0;
}
