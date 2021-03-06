/*! @file parse.cpp
 *  @brief Simple QMC Correlation function data parsing utility
 *  @author Sameed Pervaiz (pervaiz.8@osu.edu)
 *  @date 2021-05-22
 *  @copyright GPLv3
 */

// C-inherited Standard Library Headers
#include <cmath>
#include <cstdlib>

// C++ Standard Library Headers
#include <array>
#include <filesystem>
// <format> support is not always available, so defer to libfmt (same thing anyway).
// #include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>

// Third-party Libraries
#include <fmt/core.h>

// Local headers
#include "quantity.hpp"

enum cf_index
{
    MIMJ = 0,
    NINJ = 1,
    NIDJ = 2,
    DIDJ = 3,
    ZSPIN = 4,
    NUM_CFS
};

enum cf_transform
{
    TRANS_NONE = 0,
    TRANS_RECONCILE_POTENTIALS = 1
};

enum g2_type
{
    G2_NORMAL = 0,
    G2_ANTIMOMENT = 1,
    G2_VARIANTS
};

struct unaveraged_data
{
    // ( mu -> X ), where X is one of:
    //     - mu (identity)
    //     - <H> (w/ error)
    //     - <mz^2> (w/ error)
    //     - (r -> Cm(r)) (w/ error)

    // For each temperature, we handle a set of data files for various (possibly
    // repeating) chemical potentials. We'll combine the data from each file on
    // the fly here, and perform averaging later. Each instance of this struct
    // stores a set of data corresponding to a distinct (beta, mu) pair.

    // energy samples
    std::vector<qmc::quantity<double>> energies;
    // magnetic moment samples
    std::vector<qmc::quantity<double>> moments;
    // doublon (N_up*N_down) samples
    std::vector<qmc::quantity<double>> doublons;
    // density (N_up + N_down) samples
    std::vector<qmc::quantity<double>> densities;
    // z spin ((N_up - N_down)/2) samples
    std::vector<qmc::quantity<double>> z_spins;
    // The data files provide the correlated variables we can use to calculate
    // any of our correlation functions of interest.

    // To do this, we need to take the averages of these correlations for a
    // fixed temperature, chemical potential, and radius. Afterwards, we add some
    // quadratic form in the averaged scalar variables to obtain the desired cf.

    // List of correlators:
    // ======================================
    // - moment-moment correlator <m_i^2 * m_j^2>
    // - density-density correlator <n_i * n_j>
    // - density-doublon correlator <n_i * d_j>
    //     - because of the averaging, this is equal to <n_j * d_i>
    // - doublon-doublon correlator <d_i * d_j>
    // - z-spin correlation <S_i * S_j>
    std::array<std::map<double, std::vector<qmc::quantity<double>>>, NUM_CFS> cfs;
};

struct averaged_data
{
    // It's 2:47 AM and I don't know what I'm writing anymore; there is probably
    // a better way to structure this I will think of later
    // Average energy
    qmc::quantity<double> energy;
    // Average magnetic moment
    qmc::quantity<double> moment;
    // Average doublon occupancy
    qmc::quantity<double> doublon;
    // Average density
    qmc::quantity<double> density;
    // Average z spin
    qmc::quantity<double> z_spin;
    // List of correlators:
    // averaged moment-moment correlator
    // averaged density-density correlator
    // averaged density-doublon correlator
    // averaged doublon-doublon correlator
    std::array<std::map<double, qmc::quantity<double>>, NUM_CFS> cfs;
    std::array<std::map<double, qmc::quantity<double>>, G2_VARIANTS> g2;
};

using unaveraged_map = std::map<double, unaveraged_data>;
using averaged_map = std::map<double, averaged_data>;

std::string skip_lines(std::ifstream &ofs, std::string target);
void parse_file(std::string in_fname, unaveraged_map &map);
averaged_map average_data(const unaveraged_map &map);
void print_core_data_to_file(std::string out_fname, const averaged_map &map);
void print_cf_data_to_file(std::string out_fname, const averaged_map &map, cf_index cf_idx, cf_transform xform);
void print_g2_data_to_file(std::string out_fname, const averaged_map &map, g2_type g2_variant);
void print_gdh_data_to_file(std::string out_fname, const averaged_map &map);

int main(int argc, char* argv[])
{
    if (argc < 2)
        throw std::invalid_argument("Must provide data directory!");
    namespace fs = std::filesystem;
    for (auto &p : fs::directory_iterator(argv[1]))
    {
        // The directory structure looks like U[Energy]/[Size]x[Size]/beta[Inverse Temp]/...
        // So first, let's iterate over the (inverse) temperatures
        // This looks like betaA.BC, so we just take a substring starting at char 4.
        const auto float_loc = 4u;
        auto inverse_temp = p.path().filename().string().substr(float_loc);
        /*
        if (std::stod(inverse_temp) != 1.0)
            continue;
            */
        std::cerr << "Parsing inverse temperature \u03B2 = " << inverse_temp << "\n";
        // Now, move to the specified temperature subdir so we can start parsing the
        // data files for various chemical potentials mu_i:
        fs::current_path(p);
        unaveraged_map map;
        for (auto &p : fs::directory_iterator(fs::current_path()))
        {
            // For now, only examine files from the first realization (*r1.out),
            // and only files that contain the local moment correlations (r*.out)
            // TODO: add support for multiple realizations, and average over them.
            auto str = p.path().filename().string();
            if (str.at(0) == 'r')
            {
                parse_file(str, map);
            }
        }
        const auto final_temp_map = average_data(map);
        fs::current_path(fs::current_path().parent_path().parent_path());
        print_core_data_to_file("beta_" + inverse_temp + "_core.dat", final_temp_map);
        print_cf_data_to_file  ("beta_" + inverse_temp + "_cm.dat",   final_temp_map, MIMJ, TRANS_NONE);
        print_cf_data_to_file  ("beta_" + inverse_temp + "_cs.dat",   final_temp_map, ZSPIN, TRANS_NONE);
        print_cf_data_to_file  ("beta_" + inverse_temp + "_nn.dat",   final_temp_map, NINJ, TRANS_NONE);
        print_cf_data_to_file  ("beta_" + inverse_temp + "_nd.dat",   final_temp_map, NIDJ, TRANS_NONE);
        print_cf_data_to_file  ("beta_" + inverse_temp + "_dd.dat",   final_temp_map, DIDJ, TRANS_NONE);
        print_g2_data_to_file  ("beta_" + inverse_temp + "_g2.dat",   final_temp_map, G2_NORMAL);
        print_g2_data_to_file  ("beta_" + inverse_temp + "_g2bar.dat",final_temp_map, G2_ANTIMOMENT);
        print_gdh_data_to_file ("beta_" + inverse_temp + "_gdh.dat",  final_temp_map);
    }
    return 0;
}

std::string skip_lines(std::ifstream &ofs, std::string target)
{
    std::string line;
    do { std::getline(ofs, line); } while (line.find(target) == line.npos);
    return line;
}

void parse_file(std::string in_fname, unaveraged_map &map)
{
    // Now, extract the chemical potential, I guess
    std::stringstream sstream(in_fname.substr(in_fname.find("mu")+2));
    double mu;
    sstream >> mu;
    if (in_fname.find("nr") != in_fname.npos)
    {
        mu = -std::fabs(mu);
    }
    std::cerr << "Parsing chemical potential \u03BC = " << mu << std::endl;

    std::ifstream data_file(in_fname);

    // Now, parse the file.
    std::string line;

    // First, n_up
    line = skip_lines(data_file, "Average density");

    sstream = std::stringstream(line.substr(line.find("=")+2));
    qmc::quantity<double> density;
    sstream >> density.value() >> density.uncertainty();

    map[mu].densities.push_back(density);

    // Now, up occupancy and down occupancy, which we use for the average spin
    // By linearity, <Sz> = <(nup-ndown)/2>
    //                    = <nup - ndown>/2
    //                    = (<nup> - <ndown>)/2

    qmc::quantity<double> nup, ndown;

    std::getline(data_file, line);
    sstream = std::stringstream(line.substr(line.find("=")+2));
    sstream >> nup.value() >> nup.uncertainty();

    std::getline(data_file, line);
    sstream = std::stringstream(line.substr(line.find("=")+2));
    sstream >> ndown.value() >> ndown.uncertainty();

    map[mu].z_spins.push_back((nup - ndown)/2);

    // Now, the energy
    line = skip_lines(data_file, "Average Energy");
    // std::cerr << line << std::endl;

    sstream = std::stringstream(line.substr(line.find("=")+2));
    qmc::quantity<double> avg_energy;
    sstream >> avg_energy.value() >> avg_energy.uncertainty();
    // Sometimes, energies are **********. so just ignore them for now?
    if (sstream.fail())
    {
        std::cerr << "NaN energy found!\n";
        avg_energy.value() = std::nan("");
        avg_energy.uncertainty() = std::nan("");
        sstream.clear();
    }
    map[mu].energies.push_back(avg_energy);

    // Doublon number
    line = skip_lines(data_file, "Average Nup*Ndn");

    // std::cerr << "Average d line: " << line << std::endl;
    sstream = std::stringstream(line.substr(line.find("=")+2));
    qmc::quantity<double> avg_nup_ndown;
    sstream >> avg_nup_ndown.value() >> avg_nup_ndown.uncertainty();

    // std::cerr << "Found average doublon occupation: " << avg_nup_ndown.value() << std::endl;
    map[mu].doublons.push_back(avg_nup_ndown);

    // std::cerr << "Density-doublon uncorrelated part: " << (density * avg_nup_ndown).value() << std::endl;
    // std::cerr << "Skipping non-density-density correlation data!\n";
    skip_lines(data_file, "density-density correlation fn");

    // Multiple distinct triples may have the same hypotenuse length,
    // even among pythagorean triples (which are all integers). We have
    // to average over these, so we must store a vector of correlations
    // and errors for each distance, and perform the (weighted) averaging
    // later.

    // Skip the header line and start reading the cf data into our map
    std::getline(data_file, line);
    while (line.find("nud-nud") == line.npos)
    {
        // std::cerr << "Density-density line found: " << line << std::endl;
        double x, y;
        qmc::quantity<double> nup_nup, nup_ndown;
        sstream = std::stringstream(line);
        sstream >> x >> y;
        const auto r = std::hypot(x, y);

        sstream >> nup_nup.value();
        sstream.ignore(std::numeric_limits<std::streamsize>::max(), '-');
        sstream >> nup_nup.uncertainty();

        sstream >> nup_ndown.value();
        sstream.ignore(std::numeric_limits<std::streamsize>::max(), '-');
        sstream >> nup_ndown.uncertainty();

        // Because of the symmetry between up spins and down spins, the up-up
        // correlations should be equal to the down-down correlations. But I'm
        // not entirely sure of this, so let's double-check this later after
        // initial plots.

        const auto ninj = 2. * (nup_nup + nup_ndown);

        map[mu].cfs[NINJ][r].push_back(ninj);

        std::getline(data_file, line);
    }

    // Immediately after the density-density cf is the doublon-doublon density cf
    std::getline(data_file, line);
    while (line.find("0 0 dd") == line.npos)
    {
        // std::cerr << "Doublon-doublon line found: " << line << std::endl;
        sstream = std::stringstream(line);
        double x, y;
        qmc::quantity<double> didj;

        sstream >> x >> y >> didj.value() >> didj.uncertainty();
        const auto r = std::hypot(x, y);

        map[mu].cfs[DIDJ][r].push_back(didj);
        std::getline(data_file, line);
    }

    skip_lines(data_file, "mi2x-mi2x correlation function");

    // Skip the header line again
    std::getline(data_file, line);
    while (line.find("local") == line.npos)
    {
        // std::cerr << "Moment-moment line found: " << line << std::endl;
        double x, y;
        qmc::quantity<double> mimj;
        sstream = std::stringstream(line);
        sstream >> x >> y >> mimj.value() >> mimj.uncertainty();
        const double r = std::hypot(x, y);
        // std::cerr << "Line corresponds to radius " << r << std::endl;
        // the entries are coordinate symmetric, so ignore the lower
        // triangle of the (x,y) matrix
        if (y > x)
        {
            // TODO: organize the logic better
            // std::cerr << "y > x, so skipping:\n";
            std::getline(data_file, line);
            continue;
        }

        // Now, tack on our values at the given radius
        map[mu].cfs[MIMJ][r].push_back(mimj);

        std::getline(data_file, line);
    }

    // Now, density-doublon correlation
    skip_lines(data_file, "nup-nud correlation function:");
    std::getline(data_file, line);

    std::vector<std::vector<qmc::quantity<double>>> n_nuds;
    // First, put in the nup-nud values
    // TODO: do this properly, kind of just a stop-gap solution for now
    int last_index = -1;
    while (line.find("0 0 nup.d") == line.npos)
    {
        // std::cerr << "up-doublon line found: " << line << std::endl;
        sstream = std::stringstream(line);
        int x, y;
        qmc::quantity<double> nup_nud;
        sstream >> x >> y;
        sstream >> nup_nud.value() >> nup_nud.uncertainty();
        if (x > last_index)
        {
            last_index++;
            n_nuds.push_back(std::vector<qmc::quantity<double>>());
        }
        n_nuds.at(x).push_back(nup_nud);

        // std::cerr << "Pushed back values " << n_nuds.at(x).at(y).value() << " "
        //    << n_nuds.at(x).at(y).uncertainty() << std::endl;

        std::getline(data_file, line);
    }

    // Now, the other half of nidjs
    skip_lines(data_file, "nd-nud correlation function:");
    std::getline(data_file, line);

    while (line.find("0 0 ndn.d") == line.npos)
    {
        // std::cerr << "down-doublon line found: " << line << std::endl;
        sstream = std::stringstream(line);
        int x, y;
        qmc::quantity<double> nd_nud;
        sstream >> x >> y >> nd_nud.value() >> nd_nud.uncertainty();

        n_nuds.at(x).at(y) += nd_nud;
        // std::cerr << "Found sum values " << n_nuds.at(x).at(y).value() << " "
        //    << n_nuds.at(x).at(y).uncertainty() << std::endl;
        std::getline(data_file, line);
    }

    // Now, move it to our map
    for (auto x = 0u; x < n_nuds.size(); x++)
    {
        for (auto y = 0u; y < n_nuds.at(x).size(); y++)
        {
            const auto r = std::hypot(x, y);
            // std::cerr << "Pushing correlation value " << n_nuds.at(x).at(y).value() << " at radius " << r << std::endl;
            map[mu].cfs[NIDJ][r].push_back(n_nuds.at(x).at(y));
        }
    }

    // Now, the zz spin cf
    skip_lines(data_file, "zz Spin correlation");
    std::getline(data_file, line);
    while (line.find("local moment zz=") == line.npos)
    {
        double x, y;
        qmc::quantity<double> z_spin;
        sstream = std::stringstream(line);
        sstream >> x >> y >> z_spin.value() >> z_spin.uncertainty();
        const auto r = std::hypot(x,y);
        if (std::isnan(z_spin.value()) || std::isnan(z_spin.uncertainty()))
            std::exit(12);
        if (y >= x)
        {
            // std::cerr << r << std::endl;
            map[mu].cfs[ZSPIN][r].push_back(z_spin);
        }
        std::getline(data_file, line);
    }

    // Now, local moment.
    sstream = std::stringstream(line.substr(line.find("=")+1));
    qmc::quantity<double> local_moment;
    sstream >> local_moment.value() >> local_moment.uncertainty();
    // std::cerr << "Chemical potential " << mu << " has local moment " << local_moment << std::endl;

    // We will square AFTER we average all the local moments over the number
    // of realizations. This might be less numerically desirable than the original
    // ordering, but let's revisit this later.

    map[mu].moments.push_back(local_moment);
}

averaged_map average_data(const unaveraged_map &map)
{
    averaged_map final_map;

    // Okay, so for a given temperature and chemical potential, we have basically
    // turned a sequence of files with one entry per variable into a single file
    // with a sequence of data entries per variable. Now we have to average over
    // those files (in other words, over these new data entries).

    // So first, iterate over these chemical potentials, since the chemical potential
    // files are stored in the same directory (if they were in separate directories,
    // our lives would be easier)

    // Loop over chemical potential map
    for (auto it = map.begin(); it != map.end(); it++)
    {
        const auto mu = it->first;
        const auto &data = it->second;
        auto &mudata = final_map[mu];
        // Now, accumulate energies and moments, which have the same length.
        // The energies, moments, etc. are weighted averages of the respeccted
        // vectors, and the errors are handled by qmc::quantity.

        mudata.energy = qmc::average(data.energies);
        mudata.moment = qmc::average(data.moments);
        mudata.doublon = qmc::average(data.doublons);
        mudata.density = qmc::average(data.densities);
        mudata.z_spin = qmc::average(data.z_spins);

        // Now, for each radius, we perform averaging of the corr. fun, and
        // add the errors in quadrature

        for (auto i = 0; i < NUM_CFS; i++)
        {
            for (auto jt = data.cfs.at(i).begin(); jt != data.cfs.at(i).end(); jt++)
            {
                const auto r = jt->first;
                const auto &cf = jt->second;
                mudata.cfs[i][r] = qmc::average(cf);
            }
        }

        // Now, use the averages to calculate g2 and g2bar
        for (auto it = mudata.cfs.at(MIMJ).begin(); it != mudata.cfs.at(MIMJ).end(); it++)
        {
            const auto r = it->first;
            const auto g2 = it->second / (mudata.moment * mudata.moment);
            const auto g2bar = (1. - 2. * mudata.moment + it->second) /
                            ((1. - mudata.moment) * (1. - mudata.moment));
            mudata.g2[G2_NORMAL][r] = g2;
            mudata.g2[G2_ANTIMOMENT][r] = g2bar;
        }
    }
    return final_map;
}

void print_core_data_to_file(std::string out_fname, const averaged_map &map)
{
    std::ofstream output_file(out_fname);
    // Now, print to files
    output_file << "# mu   <H>    delta_<H>    <mz^2>    delta_<mz^2>    <n>    delta_<n>    <d>    delta_<d>    <nd>    delta_<nd>    ";
    // Take the first chemical potential, and print out all the lattice distances
    // it doesn't matter which potential we use, since they're all on the same lattice.
    output_file << std::endl;

    // Print the core data as a fn of mu
    for (auto it = map.begin(); it != map.end(); it++)
    {
        const auto mu = it->first;
        const auto &data = it->second;
        output_file << mu << " " << data.energy.value() << " " << data.energy.uncertainty() << " "
            << data.moment.value() << " " << data.moment.uncertainty() << " "
            << data.density.value() << " " << data.density.uncertainty() << " "
            << data.doublon.value() << " " << data.doublon.uncertainty();

        output_file << std::endl;
    }
}

void print_cf_data_to_file(std::string out_fname, const averaged_map &map, cf_index cf_idx, cf_transform xform)
{
    // std::ofstream debug_out("debug_" + out_fname);
    std::ofstream output_file(out_fname);
    // Now, print to files
    output_file << "# mu   <H>    delta_<H>    <mz^2>    delta_<mz^2>    <n>    delta_<n>    ";
    for (auto it = map.begin()->second.cfs.at(cf_idx).begin(); it != map.begin()->second.cfs.at(cf_idx).end(); it++)
    {
        const auto r = it->first;
        output_file << "CF(" << std::setprecision(3) << r <<
            ")    dCF(" << std::setprecision(3) << r << ")    ";
    }
    output_file << std::endl;
    output_file << std::setprecision(4);

    // Print the cf data as a fn of mu
    for (auto it = map.begin(); it != map.end(); it++)
    {
        const auto mu = it->first;
        auto &data = it->second;
        const std::array<qmc::quantity<double>, NUM_CFS> quad_vals
        ({
            data.moment * data.moment,
            data.density * data.density,
            data.density * data.doublon,
            data.doublon * data.doublon,
            data.z_spin * data.z_spin
        });

        output_file << mu << " " << data.energy.value() << " " << data.energy.uncertainty() << " "
            << data.moment.value() << " " << data.moment.uncertainty() << " "
            << data.density.value() << " " << data.density.uncertainty() << " ";

        for (auto jt = data.cfs.at(cf_idx).begin(); jt != data.cfs.at(cf_idx).end(); jt++)
        {
            auto r = jt->first;
            // If we're above half-filling, we need to invert the sites to obtain
            // data comparable to <n> less than 1. This leads to different expressions
            // for the correlation functions, which have been calculated elsewhere.

            // However, exceptions exist, like for Cm, which is invariant under this
            // inversion, despite its components not being such.

            // Maybe we could consider including the reasoning here, too.
            auto cf = jt->second - quad_vals.at(cf_idx);

            /*
            if (r == 1)
                debug_out << std::setw(12) << mu
                    << std::setw(12) << jt->second
                    << std::setw(12) << quad_vals.at(cf_idx)
                    << std::setw(12) << cf << std::endl;
                    */
            if (xform == TRANS_RECONCILE_POTENTIALS && mu > 0)
            {
                if (cf_idx == NIDJ)
                {
                    // <n'd'> = <nn> - <nd>
                    cf = data.cfs.at(NINJ).at(r) - quad_vals.at(NINJ) - cf;
                }
                else if (cf_idx == DIDJ)
                {
                    cf += data.cfs.at(NINJ).at(r) - quad_vals.at(NINJ)
                        - 2. * (data.cfs.at(NIDJ).at(r) - quad_vals.at(NIDJ));
                }
            }
            output_file << cf << " ";
        }
        output_file << std::endl;
    }
}

void print_g2_data_to_file(std::string out_fname, const averaged_map &map, g2_type g2_variant)
{
    std::ofstream output_file(out_fname);
    output_file << "# mu   <H>    delta_<H>    <mz^2>    delta_<mz^2>    <n>    delta_<n>    ";
    for (auto it = map.begin()->second.g2.at(g2_variant).begin(); it != map.begin()->second.g2.at(g2_variant).end(); it++)
    {
        const auto r = it->first;
        output_file << "g2(" << std::setprecision(3) << r <<
            ")    dg2(" << std::setprecision(3) << r << ")    ";
    }
    output_file << std::endl;
    output_file << std::setprecision(4);

    for (auto it = map.begin(); it != map.end(); it++)
    {
        const auto mu = it->first;
        auto &data = it->second;

        output_file << mu << " " << data.energy.value() << " " << data.energy.uncertainty() << " "
            << data.moment.value() << " " << data.moment.uncertainty() << " "
            << data.density.value() << " " << data.density.uncertainty() << " ";

        for (auto jt = data.g2.at(g2_variant).begin(); jt != data.g2.at(g2_variant).end(); jt++)
        {
            output_file << jt->second << " ";
        }
        output_file << std::endl;
    }
}

void print_gdh_data_to_file(std::string out_fname, const averaged_map &map)
{
    std::ofstream output_file(out_fname);
    output_file << "# mu   <H>    delta_<H>    <mz^2>    delta_<mz^2>    <n>    delta_<n>    ";
    for (auto it = map.begin()->second.cfs.at(NIDJ).begin(); it != map.begin()->second.cfs.at(NIDJ).end(); it++)
    {
        const auto r = it->first;
        output_file << "gdh(" << std::setprecision(3) << r <<
            ")    dgh(" << std::setprecision(3) << r << ")    ";
    }
    output_file << std::endl;
    output_file << std::setprecision(4);

    for (auto it = map.begin(); it != map.end(); it++)
    {
        const auto mu = it->first;
        auto &data = it->second;

        output_file << mu << " " << data.energy.value() << " " << data.energy.uncertainty() << " "
            << data.moment.value() << " " << data.moment.uncertainty() << " "
            << data.density.value() << " " << data.density.uncertainty() << " ";

        for (auto jt = data.cfs.at(NIDJ).begin(); jt != data.cfs.at(NIDJ).end(); jt++)
        {
            // <dh> =  <d(2-n)>
            //      = 2<d>-<dn>
            output_file << (2.*data.doublon - jt->second)/(data.doublon * (2.-data.density)) << " ";
        }
        output_file << std::endl;
    }
}
