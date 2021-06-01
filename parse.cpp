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

// TODO: consider map of structs, rather than struct of maps
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
    std::vector<double> energies;
    std::vector<double> d_energies;
    // magnetic moment samples
    std::vector<double> moments;
    std::vector<double> d_moments;
    // doublon (N_up*N_down) samples
    std::vector<double> doublons;
    std::vector<double> d_doublons;
    // density (N_up + N_down) samples
    std::vector<double> densities;
    std::vector<double> d_densities;
    // The data files provide the correlated variables we can use to calculate
    // any of our correlation functions of interest.

    // To do this, we need to take the averages of these correlations for a
    // fixed temperature, chemical potential, and radius. Afterwards, we add some
    // quadratic form in the averaged scalar variables to obtain the desired cf.
    // moment-moment correlator <m_i * m_j>
    std::map<double, std::vector<double>> mimjs;
    std::map<double, std::vector<double>> d_mimjs;
    // density-density correlator <n_i * n_j>
    std::map<double, std::vector<double>> ninjs;
    std::map<double, std::vector<double>> d_ninjs;
    // density-doublon correlator <n_i * d_j>
    // because of the averaging, this is equal to <n_j * d_i>
    std::map<double, std::vector<double>> nidjs;
    std::map<double, std::vector<double>> d_nidjs;
    // doublon-doublon correlator <d_i * d_j>
    std::map<double, std::vector<double>> didjs;
    std::map<double, std::vector<double>> d_didjs;
};

struct averaged_data
{
    // It's 2:47 AM and I don't know what I'm writing anymore; there is probably
    // a better way to structure this I will think of later
    // Average energy
    double energy;
    double d_energy;
    // Average magnetic moment
    double moment;
    double d_moment;
    // Average magnetic moment squared
    double moment_sq;
    double d_moment_sq;
    // Average doublon occupancy
    double doublon;
    double d_doublon;
    // Average doublon occupancy squared
    double doublon_sq;
    double d_doublon_sq;
    // Average density
    double density;
    double d_density;
    // Average density squared
    double density_sq;
    double d_density_sq;
    // Average doublon occupancy * density
    double density_doublon;
    double d_density_doublon;
    // averaged moment-moment correlator
    std::map<double, double> mimj;
    std::map<double, double> d_mimj;
    // averaged density-density correlator
    std::map<double, double> ninj;
    std::map<double, double> d_ninj;
    // averaged density-doublon correlator
    std::map<double, double> nidj;
    std::map<double, double> d_nidj;
    // averaged doublon-doublon correlator
    std::map<double, double> didj;
    std::map<double, double> d_didj;
};

using unaveraged_map = std::map<double, unaveraged_data>;
using averaged_map = std::map<double, averaged_data>;

void parse_file(std::string in_fname, unaveraged_map &map);
averaged_map average_data(const unaveraged_map &map);
void print_core_data_to_file(std::string out_fname, const averaged_map &map);
void print_cm_data_to_file(std::string out_fname, const averaged_map &map);
void print_dd_data_to_file(std::string out_fname, const averaged_map &map);
void print_dh_data_to_file(std::string out_fname, const averaged_map &map);
void print_hh_data_to_file(std::string out_fname, const averaged_map &map);

int main()
{
    namespace fs = std::filesystem;
    for (auto &p : fs::directory_iterator("10x10"))
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
        print_cm_data_to_file  ("beta_" + inverse_temp + "_cm.dat",   final_temp_map);
        print_dd_data_to_file  ("beta_" + inverse_temp + "_dd.dat",   final_temp_map);
        print_dh_data_to_file  ("beta_" + inverse_temp + "_dh.dat",   final_temp_map);
        print_hh_data_to_file  ("beta_" + inverse_temp + "_hh.dat",   final_temp_map);
    }
    return 0;
}

void parse_file(std::string in_fname, unaveraged_map &map)
{
    // Now, extract the chemical potential, I guess
    std::stringstream sstream(in_fname.substr(in_fname.find("mu")+2));
    double mu;
    sstream >> mu;
    /*
    if (mu != 10.)
        return;
        */
    std::cerr << "Parsing chemical potential \u03BC = " << mu << std::endl;

    std::ifstream data_file(in_fname);

    // Now, parse the file.
    std::string line;

    // First, n_up
    do
    {
        std::getline(data_file, line);
    }
    while (line.find("Average density") == line.npos);

    sstream = std::stringstream(line.substr(line.find("=")+2));
    double density, density_uc;
    sstream >> density >> density_uc;

    map[mu].densities.push_back(density);
    map[mu].d_densities.push_back(density_uc);

    // Now, the energy
    do
    {
        std::getline(data_file, line);
    }
    while (line.find("Average Energy") == line.npos);
    std::cerr << line << std::endl;

    sstream = std::stringstream(line.substr(line.find("=")+2));
    double avg_energy, energy_uc;
    sstream >> avg_energy >> energy_uc;
    // Sometimes, energies are **********. so just ignore them for now?
    if (sstream.fail())
    {
        std::cerr << "NaN energy found!\n";
        avg_energy = std::nan("");
        energy_uc = std::nan("");
        sstream.clear();
    }
    map[mu].energies.push_back(avg_energy);
    map[mu].d_energies.push_back(energy_uc);

    // Doublon number
    do
    {
        std::getline(data_file, line);
    }
    while (line.find("Average Nup*Ndn") == line.npos);

    std::cerr << "Average d line: " << line << std::endl;
    sstream = std::stringstream(line.substr(line.find("=")+2));
    double avg_nup_ndown, nup_ndown_uc;
    sstream >> avg_nup_ndown >> nup_ndown_uc;

    std::cerr << "Found average doublon occupation: " << avg_nup_ndown << std::endl;
    map[mu].doublons.push_back(avg_nup_ndown);
    map[mu].d_doublons.push_back(nup_ndown_uc);

    std::cerr << "Density-doublon uncorrelated part: " << density * avg_nup_ndown << std::endl;
    std::cerr << "Skipping non-density-density correlation data!\n";
    do
    {
        std::getline(data_file, line);
    }
    while (line.find("density-density correlation fn") == line.npos);

    // Multiple distinct triples may have the same hypotenuse length,
    // even among pythagorean triples (which are all integers). We have
    // to average over these, so we must store a vector of correlations
    // and errors for each distance, and perform the (weighted) averaging
    // later.

    // Skip the header line and start reading the cf data into our map
    std::getline(data_file, line);
    while (line.find("nud-nud") == line.npos)
    {
        std::cerr << "Density-density line found: " << line << std::endl;
        double x, y,
               nup_nup, nup_nup_uc,
               nup_ndown, nup_ndown_uc;
        sstream = std::stringstream(line);
        sstream >> x >> y;
        const auto r = std::hypot(x, y);

        sstream >> nup_nup;
        sstream.ignore(std::numeric_limits<std::streamsize>::max(), '-');
        sstream >> nup_nup_uc;

        sstream >> nup_ndown;
        sstream.ignore(std::numeric_limits<std::streamsize>::max(), '-');
        sstream >> nup_ndown_uc;

        // Because of the symmetry between up spins and down spins, the up-up
        // correlations should be equal to the down-down correlations. But I'm
        // not entirely sure of this, so let's double-check this later after
        // initial plots.

        const auto ninj = 2. * (nup_nup + nup_ndown);
        const auto ninj_uc = 2. * std::sqrt(nup_nup_uc * nup_nup_uc +
                                            nup_ndown_uc * nup_ndown_uc);

        map[mu].ninjs[r].push_back(ninj);
        map[mu].d_ninjs[r].push_back(ninj_uc);

        std::getline(data_file, line);
    }

    // Immediately after the density-density cf is the doublon-doublon density cf
    std::getline(data_file, line);
    while (line.find("0 0 dd") == line.npos)
    {
        std::cerr << "Doublon-doublon line found: " << line << std::endl;
        sstream = std::stringstream(line);
        double x, y,
               didj, didj_uc;

        sstream >> x >> y >> didj >> didj_uc;
        const auto r = std::hypot(x, y);

        map[mu].didjs[r].push_back(didj);
        map[mu].d_didjs[r].push_back(didj_uc);
        std::getline(data_file, line);
    }

    do
    {
        std::getline(data_file, line);
    }
    while (line.find("mi2x-mi2x correlation function") == line.npos);

    // Skip the header line again
    std::getline(data_file, line);
    while (line.find("local") == line.npos)
    {
        std::cerr << "Moment-moment line found: " << line << std::endl;
        double x, y,
               mimj, mimj_uc;
        sstream = std::stringstream(line);
        sstream >> x >> y >> mimj >> mimj_uc;
        const double r = std::hypot(x, y);
        std::cerr << "Line corresponds to radius " << r << std::endl;
        // the entries are coordinate symmetric, so ignore the lower
        // triangle of the (x,y) matrix
        if (y > x)
        {
            // TODO: organize the logic better
            std::cerr << "y > x, so skipping:\n";
            std::getline(data_file, line);
            continue;
        }

        // Now, tack on our values at the given radius
        map[mu].mimjs[r].push_back(mimj);
        map[mu].d_mimjs[r].push_back(mimj_uc);

        std::getline(data_file, line);
    }

    // Now, density-doublon correlation
    do
    {
        std::getline(data_file, line);
    }
    while (line.find("nup-nud correlation function:") == line.npos);
    std::getline(data_file, line);

    std::vector<std::vector<double>> n_nuds, n_nuds_uc;
    // First, put in the nup-nud values
    // TODO: do this properly, kind of just a stop-gap solution for now
    int last_index = -1;
    while (line.find("0 0 nup.d") == line.npos)
    {
        std::cerr << "up-doublon line found: " << line << std::endl;
        sstream = std::stringstream(line);
        int x, y;
        double nup_nud, nup_nud_uc;
        sstream >> x >> y;
        sstream >> nup_nud >> nup_nud_uc;
        if (x > last_index)
        {
            last_index++;
            n_nuds.push_back(std::vector<double>());
            n_nuds_uc.push_back(std::vector<double>());
        }
        n_nuds.at(x).push_back(nup_nud);
        n_nuds_uc.at(x).push_back(nup_nud_uc);

        std::cerr << "Pushed back values " << n_nuds.at(x).at(y) << " " << n_nuds_uc.at(x).at(y) << std::endl;

        std::getline(data_file, line);
    }

    // Now, the other half of nidjs
    do
    {
        std::getline(data_file, line);
    }
    while (line.find("nd-nud correlation function:") == line.npos);
    std::getline(data_file, line);

    while (line.find("0 0 ndn.d") == line.npos)
    {
        std::cerr << "down-doublon line found: " << line << std::endl;
        sstream = std::stringstream(line);
        int x, y;
        double nd_nud, nd_nud_uc;
        sstream >> x >> y >> nd_nud >> nd_nud_uc;

        n_nuds.at(x).at(y) += nd_nud;
        n_nuds_uc.at(x).at(y) = std::sqrt(n_nuds_uc.at(x).at(y) * n_nuds_uc.at(x).at(y) + nd_nud_uc * nd_nud_uc);
        std::cerr << "Found sum values " << n_nuds.at(x).at(y) << " " << n_nuds_uc.at(x).at(y) << std::endl;
        std::getline(data_file, line);
    }

    // Now, move it to our map
    for (auto x = 0u; x < n_nuds.size(); x++)
    {
        for (auto y = 0u; y < n_nuds.at(x).size(); y++)
        {
            const auto r = std::hypot(x, y);
            std::cout << "Pushing correlation value " << n_nuds.at(x).at(y) << " at radius " << r << std::endl;
            map[mu].nidjs[r].push_back(n_nuds.at(x).at(y));
            map[mu].d_nidjs[r].push_back(n_nuds_uc.at(x).at(y));
        }
    }

    // Now, local moment.
    do
    {
        std::getline(data_file, line);
    }
    while (line.find("0 0 local moment zz=") == line.npos);

    sstream = std::stringstream(line.substr(line.find("=")+1));
    double local_moment, local_moment_uc;
    sstream >> local_moment >> local_moment_uc;

    // We will square AFTER we average all the local moments over the number
    // of realizations. This might be less numerically desirable than the original
    // ordering, but let's revisit this later.

    map[mu].moments.push_back(local_moment);
    map[mu].d_moments.push_back(local_moment_uc);
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
        // The energies and moments are simply the averages of the corresponding
        // vectors, so we use std::accumulate to take advantage of parallelism (note:
        // this might require explicitly detailing the execution policy)

        // NOTE: availability of std::reduce and std::transform_reduce is absolutely
        // awful as of now (26-May-2021), so we'll stick to the serial versions
        // for now...

        // For the errors, we add them in quadrature and then divide by N. We
        // can accomplish this quickly by taking the inner product of the vector
        // with itself, and then manipulating the result.
        mudata.energy = std::accumulate(data.energies.begin(), data.energies.end(), 0.) / data.energies.size();
        mudata.d_energy =
            std::sqrt(std::inner_product(data.d_energies.begin(), data.d_energies.end(), data.d_energies.begin(), 0.)) / data.energies.size();

        mudata.moment = std::accumulate(data.moments.begin(), data.moments.end(), 0.) / data.moments.size();
        mudata.d_moment =
            std::sqrt(std::inner_product(data.d_moments.begin(), data.d_moments.end(), data.d_moments.begin(), 0.)) / data.moments.size();

        mudata.doublon = std::accumulate(data.doublons.begin(), data.doublons.end(), 0.) / data.doublons.size();
        mudata.d_doublon =
            std::sqrt(std::inner_product(data.d_doublons.begin(), data.d_doublons.end(), data.d_doublons.begin(), 0.)) / data.doublons.size();

        mudata.density = std::accumulate(data.densities.begin(), data.densities.end(), 0.) / data.densities.size();
        mudata.d_density =
            std::sqrt(std::inner_product(data.d_densities.begin(), data.d_densities.end(), data.d_densities.begin(), 0.)) / data.densities.size();

        mudata.density_doublon = mudata.doublon * mudata.density;
        mudata.d_density_doublon = mudata.density_doublon *
            std::sqrt((mudata.d_density / mudata.density) * (mudata.d_density / mudata.density) +
                    (mudata.d_doublon / mudata.doublon) * (mudata.d_doublon / mudata.doublon));
        // So far we've only looked at the local moment, but we want the squared
        // local moment, so we recalculate the value and its error

        // For the uncertainty, recall the following:
        // delta(a*b) = |a*b| sqrt( (delta(a)/a)^2 + (delta(b)/b)^2 )
        // In the case that a=b, this simplifies to
        // delta(a^2) = |a^2| sqrt( 2 (delta(a)/a)^2 )
        //            = |a^2| |delta(a)/a| sqrt(2)
        //            = |a * delta(a)| sqrt(2)
        mudata.moment_sq = mudata.moment * mudata.moment;
        mudata.d_moment_sq = std::abs(mudata.moment * mudata.d_moment) * std::sqrt(2.);

        mudata.density_sq = mudata.density * mudata.density;
        mudata.d_density_sq = std::abs(mudata.density * mudata.d_density) * std::sqrt(2.);

        mudata.doublon_sq = mudata.doublon * mudata.doublon;
        mudata.d_doublon_sq = std::abs(mudata.doublon * mudata.d_doublon) * std::sqrt(2.);

        // Now, for each radius, we perform averaging of the corr. fun, and
        // add the errors in quadrature

        for (auto jt = data.mimjs.begin(); jt != data.mimjs.end(); jt++)
        {
            const auto r = jt->first;
            const auto &mimjs = jt->second;
            mudata.mimj[r] = std::accumulate(mimjs.begin(), mimjs.end(), 0.) / mimjs.size();
            mudata.d_mimj[r] =
                std::sqrt(std::inner_product(data.d_mimjs.at(r).begin(), data.d_mimjs.at(r).end(), data.d_mimjs.at(r).begin(), 0.)) / mimjs.size();
        }

        for (auto jt = data.ninjs.begin(); jt != data.ninjs.end(); jt++)
        {
            const auto r = jt->first;
            const auto &ninjs = jt->second;
            mudata.ninj[r] = std::accumulate(ninjs.begin(), ninjs.end(), 0.) / ninjs.size();
            mudata.d_ninj[r] =
                std::sqrt(std::inner_product(data.d_ninjs.at(r).begin(), data.d_ninjs.at(r).end(), data.d_ninjs.at(r).begin(), 0.)) / ninjs.size();
        }

        for (auto jt = data.nidjs.begin(); jt != data.nidjs.end(); jt++)
        {
            const auto r = jt->first;
            const auto &nidjs = jt->second;
            mudata.nidj[r] = std::accumulate(nidjs.begin(), nidjs.end(), 0.) / nidjs.size();
            mudata.d_nidj[r] =
                std::sqrt(std::inner_product(data.d_nidjs.at(r).begin(), data.d_nidjs.at(r).end(), data.d_nidjs.at(r).begin(), 0.)) / nidjs.size();
        }

        for (auto jt = data.didjs.begin(); jt != data.didjs.end(); jt++)
        {
            const auto r = jt->first;
            const auto &didjs = jt->second;
            mudata.didj[r] = std::accumulate(didjs.begin(), didjs.end(), 0.) / didjs.size();
            mudata.d_didj[r] =
                std::sqrt(std::inner_product(data.d_didjs.at(r).begin(), data.d_didjs.at(r).end(), data.d_didjs.at(r).begin(), 0.)) / didjs.size();
        }
    }

    return final_map;
}

void print_core_data_to_file(std::string out_fname, const averaged_map &map)
{
    std::ofstream output_file(out_fname);
    // Now, print to files
    output_file << "# mu   <H>    delta_<H>    <mz^2>    delta_<mz^2>    ";
    // Take the first chemical potential, and print out all the lattice distances
    // it doesn't matter which potential we use, since they're all on the same lattice.
    output_file << std::endl;

    // Print the core data as a fn of mu
    for (auto it = map.begin(); it != map.end(); it++)
    {
        const auto mu = it->first;
        const auto &data = it->second;
        output_file << mu << " " << data.energy << " " << data.d_energy << " "
            << data.moment_sq << " " << data.d_moment_sq << " ";

        output_file << std::endl;
    }
}

void print_cm_data_to_file(std::string out_fname, const averaged_map &map)
{
    std::ofstream output_file(out_fname);
    // Now, print to files
    output_file << "# mu   <H>    delta_<H>    <mz^2>    delta_<mz^2>    ";
    for (auto it = map.begin()->second.mimj.begin(); it != map.begin()->second.mimj.end(); it++)
    {
        const auto r = it->first;
        output_file << "Cm(" << std::setprecision(3) << r <<
            ")    dCm(" << std::setprecision(3) << r << ")    ";
    }
    output_file << std::endl;

    // Print the Cm data as a fn of mu
    for (auto it = map.begin(); it != map.end(); it++)
    {
        const auto mu = it->first;
        const auto &data = it->second;
        output_file << mu << " " << data.energy << " " << data.d_energy << " "
            << data.moment_sq << " " << data.d_moment_sq << " ";

        for (auto jt1 = data.mimj.begin(), jt2 = data.d_mimj.begin();
                jt1 != data.mimj.end();
                jt1++, jt2++)
        {
            output_file << std::setprecision(8) << jt1->second - data.moment_sq << " "
                << std::sqrt(jt2->second * jt2->second + data.d_moment_sq * data.d_moment_sq) << " ";
        }
        output_file << std::endl;
    }
}

void print_dd_data_to_file(std::string out_fname, const averaged_map &map)
{
    std::ofstream output_file(out_fname);
    std::ofstream log_file (out_fname + ".log");
    // Now, print to files
    output_file << "# mu   <H>    delta_<H>    <mz^2>    delta_<mz^2>    ";
    for (auto it = map.begin()->second.didj.begin(); it != map.begin()->second.didj.end(); it++)
    {
        const auto r = it->first;
        output_file << "Cdd(" << std::setprecision(3) << r <<
            ")    dCdd(" << std::setprecision(3) << r << ")    ";
    }
    output_file << std::endl;

    // Print the Cdd data as a fn of mu
    for (auto it = map.begin(); it != map.end(); it++)
    {
        const auto mu = it->first;
        const auto &data = it->second;
        output_file << mu << " " << data.energy << " " << data.d_energy << " "
            << data.moment_sq << " " << data.d_moment_sq << " ";

        for (auto jt1 = data.didj.begin(), jt2 = data.d_didj.begin();
                jt1 != data.didj.end();
                jt1++, jt2++)
        {
            log_file << std::setprecision(8) << "Averaged didj(" << jt1->first << "): " << jt1->second
                << "\nAveraged doublon number: " << data.doublon << std::endl;
            log_file << "Correlation: " << jt1->second - data.doublon_sq << std::endl;
            output_file << std::setprecision(8) << jt1->second - data.doublon * data.doublon << " "
                << std::sqrt(jt2->second * jt2->second + data.d_doublon_sq * data.d_doublon_sq) << " ";
        }
        output_file << std::endl;
    }
}

void print_dh_data_to_file(std::string out_fname, const averaged_map &map)
{
    std::ofstream output_file(out_fname);
    std::ofstream derived_output_file("cm_derived_" + out_fname);
    std::ofstream log_file (out_fname + ".log");
    // Now, print to files
    output_file << "# mu   <H>    delta_<H>    <mz^2>    delta_<mz^2>    ";
    derived_output_file << "# mu   <H>    delta_<H>    <mz^2>    delta_<mz^2>    ";
    for (auto it = map.begin()->second.didj.begin(); it != map.begin()->second.didj.end(); it++)
    {
        const auto r = it->first;
        output_file << "Cdh(" << std::setprecision(3) << r <<
            ")    dCdh(" << std::setprecision(3) << r << ")    ";
        derived_output_file << "Cdh(" << std::setprecision(3) << r <<
            ")    dCdh(" << std::setprecision(3) << r << ")    ";
    }
    output_file << std::endl;
    derived_output_file << std::endl;

    // Print the Cdh data as a fn of mu
    for (auto it = map.begin(); it != map.end(); it++)
    {
        const auto mu = it->first;
        const auto &data = it->second;
        log_file << std::setprecision(8) << mu << std::endl;
        output_file << mu << " " << data.energy << " " << data.d_energy << " "
            << data.moment_sq << " " << data.d_moment_sq << " ";
        derived_output_file << mu << " " << data.energy << " " << data.d_energy << " "
            << data.moment_sq << " " << data.d_moment_sq << " ";

        for (auto jt1 = data.nidj.begin(), jt2 = data.d_nidj.begin(),
                  jt3 = data.mimj.begin(), jt4 = data.d_mimj.begin(),
                  jt5 = data.ninj.begin(), jt6 = data.d_ninj.begin(),
                  jt7 = data.didj.begin(), jt8 = data.d_didj.end();
                jt1 != data.nidj.end();
                jt1++, jt2++, jt3++, jt4++, jt5++, jt6++, jt7++, jt8++)
        {
            // The doublon-holon correlation is the negative doublon-density correlation.
            log_file << std::setprecision(8) << "Averaged nidj(" << jt1->first << "): " << jt1->second
                << "\nAveraged doublon number: " << data.doublon << std::endl
                << "Averaged density: " << data.density << std::endl;
            log_file << "Correlation: " << -(jt1->second - data.density_doublon) << std::endl;
            output_file << std::setprecision(8)
                <<  -(jt1->second - data.density_doublon) << " "
                << std::sqrt(jt2->second * jt2->second + data.d_density_doublon * data.d_density_doublon) << " ";
            derived_output_file << std::setprecision(8)
                <<  (jt5->second - data.density_sq) - (jt3->second - data.moment_sq) + 4. * (jt7->second - data.doublon_sq) << " "
                << std::sqrt(jt6->second * jt6->second
                        + data.d_density_sq * data.d_density_sq
                        + jt4->second * jt4->second
                        + data.d_moment_sq * data.d_moment_sq
                        + 16. * jt8->second * jt8->second
                        + 16. * data.d_doublon_sq * data.d_doublon_sq
                        ) << " ";
        }
        output_file << std::endl;
        derived_output_file << std::endl;
    }
}

void print_hh_data_to_file(std::string out_fname, const averaged_map &map)
{
    std::ofstream output_file(out_fname);
    std::ofstream log_file (out_fname + ".log");
    // Now, print to files
    output_file << "# mu   <H>    delta_<H>    <mz^2>    delta_<mz^2>    ";
    for (auto it = map.begin()->second.didj.begin(); it != map.begin()->second.didj.end(); it++)
    {
        const auto r = it->first;
        output_file << "Chh(" << std::setprecision(3) << r <<
            ")    dChh(" << std::setprecision(3) << r << ")    ";
    }
    output_file << std::endl;

    // Print the Chh data as a fn of mu
    for (auto it = map.begin(); it != map.end(); it++)
    {
        const auto mu = it->first;
        const auto &data = it->second;
        output_file << mu << " " << data.energy << " " << data.d_energy << " "
            << data.moment_sq << " " << data.d_moment_sq << " ";

        for (auto jt1 = data.ninj.begin(), jt2 = data.d_ninj.begin(), jt3 = data.d_mimj.begin();
                jt1 != data.ninj.end();
                jt1++, jt2++, jt3++)
        {
            // The holon-holon correlation is just the density-density correlation
            log_file << std::setprecision(8) << "Averaged didj(" << jt1->first << "): " << jt1->second
                << "\nAveraged density: " << data.density << std::endl;
            log_file << "Correlation: " << jt1->second - data.density_sq << std::endl;
            output_file << std::setprecision(8)
                << jt1->second - data.density_sq << " "
                << std::sqrt(jt2->second * jt2->second + data.d_density_sq * data.d_density_sq) << " ";
        }
        output_file << std::endl;
    }
}

