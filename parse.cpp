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
    std::vector<double> energies;
    std::vector<double> d_energies;
    std::vector<double> moments;
    std::vector<double> d_moments;
    std::map<double, std::vector<double>> Cms;
    std::map<double, std::vector<double>> d_Cms;
};

struct averaged_data
{
    // It's 2:47 AM and I don't know what I'm writing anymore; there is probably
    // a better way to structure this I will think of later
    double energy;
    double d_energy;
    double moment;
    double d_moment;
    std::map<double, double> Cm;
    std::map<double, double> d_Cm;
};

using unaveraged_map = std::map<double, unaveraged_data>;
using averaged_map = std::map<double, averaged_data>;

void parse_file(std::string in_fname, unaveraged_map &map);
averaged_map average_data(unaveraged_map &map);
void print_to_file(std::string out_fname, averaged_map &map);

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
        auto final_temp_map = average_data(map);
        fs::current_path(fs::current_path().parent_path().parent_path());
        print_to_file("beta_" + inverse_temp + "_parsed.dat", final_temp_map);
    }
    return 0;
}

void parse_file(std::string in_fname, unaveraged_map &map)
{
    // Now, extract the chemical potential, I guess
    std::stringstream sstream(in_fname.substr(in_fname.find("mu")+2));
    double mu;
    sstream >> mu;
    std::cerr << "Parsing chemical potential \u03BC = " << mu << std::endl;

    std::ifstream data_file(in_fname);

    // Now, parse the file. For now, let's only track the average site energy
    // So skip everything up until that line
    std::string line;
    do
    {
        std::getline(data_file, line);
    }
    while (line.find("Average Energy") == line.npos);
    std::cerr << line << std::endl;

    // Now print the line with its uncertainty as a function of chemical
    // potential
    sstream = std::stringstream(line.substr(line.find("=")+1));
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

    std::cerr << "Skipping non-moment correlation function data!\n";
    do
    {
        std::getline(data_file, line);
    }
    while (line.find("mi2x-mi2x correlation function") == line.npos);

    // Multiple distinct triples may have the same hypotenuse length,
    // even among pythagorean triples (which are all integers). We have
    // to average over these, so we must store a vector of correlations
    // and errors for each distance, and perform the (weighted) averaging
    // later.

    // Skip the header line and start reading the cf data into our map
    std::getline(data_file, line);
    while (line.find("local") == line.npos)
    {
        std::cerr << "Moment line found: " << line << std::endl;
        double x, y, cf, cf_uc;
        sstream = std::stringstream(line);
        sstream >> x >> y >> cf >> cf_uc;
        double r = std::hypot(x, y);
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
        map[mu].Cms[r].push_back(cf);
        map[mu].d_Cms[r].push_back(cf_uc);

        std::getline(data_file, line);
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

averaged_map average_data(unaveraged_map &map)
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
        auto mu = it->first;
        auto &data = it->second;
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
        final_map[mu].energy = std::accumulate(data.energies.begin(), data.energies.end(), 0.) / data.energies.size();
        final_map[mu].d_energy =
            std::sqrt(std::inner_product(data.d_energies.begin(), data.d_energies.end(), data.d_energies.begin(), 0.)) / data.energies.size();

        final_map[mu].moment = std::accumulate(data.moments.begin(), data.moments.end(), 0.) / data.moments.size();
        final_map[mu].d_moment =
            std::sqrt(std::inner_product(data.d_moments.begin(), data.d_moments.end(), data.d_moments.begin(), 0.)) / data.moments.size();

        // So far we've only looked at the local moment, but we want the squared
        // local moment, so we recalculate the value and its error

        // For the uncertainty, recall the following:
        // delta(a*b) = |a*b| sqrt( (delta(a)/a)^2 + (delta(b)/b)^2 )
        // In the case that a=b, this simplifies to
        // delta(a^2) = |a^2| sqrt( 2 (delta(a)/a)^2 )
        //            = |a^2| |delta(a)/a| sqrt(2)
        //            = |a * delta(a)| sqrt(2)
        final_map[mu].d_moment = std::abs(final_map[mu].moment * final_map[mu].d_moment) * std::sqrt(2.);
        final_map[mu].moment *= final_map[mu].moment;


        // Now, for each radius, we perform averaging of the corr. fun, and
        // add the errors in quadrature

        for (auto jt = data.Cms.begin(); jt != data.Cms.end(); jt++)
        {
            auto r = jt->first;
            auto &cms = jt->second;
            final_map[mu].Cm[r] = std::accumulate(cms.begin(), cms.end(), 0.) / cms.size();
            final_map[mu].d_Cm[r] = std::sqrt(std::inner_product(data.d_Cms[r].begin(), data.d_Cms[r].end(), data.d_Cms[r].begin(), 0.)) / cms.size();
        }
    }

    return final_map;
}

void print_to_file(std::string out_fname, averaged_map &map)
{
    std::ofstream output_file(out_fname);
    // Now, print to files
    output_file << "# mu   <H>    delta_<H>    <mz^2>    delta_<mz^2>    ";
    // Take the first chemical potential, and print out all the lattice distances
    // it doesn't matter which potential we use, since they're all on the same lattice.
    for (auto it = map.begin()->second.Cm.begin(); it != map.begin()->second.Cm.end(); it++)
    {
        const auto r = it->first;
        output_file << "Cm(" << std::setprecision(3) << r <<
            ")    dCm(" << std::setprecision(3) << r << ")    ";
    }
    output_file << std::endl;

    // Now, loop over chemical potentials and print the data for each
    for (auto it = map.begin(); it != map.end(); it++)
    {
        const auto mu = it->first;
        const auto &data = it->second;
        output_file << mu << " " << data.energy << " " << data.d_energy << " "
            << data.moment << " " << data.d_moment << " ";

        for (auto jt1 = data.Cm.begin(), jt2 = data.d_Cm.begin();
                jt1 != data.Cm.end();
                jt1++, jt2++)
        {
            output_file << std::setprecision(8) << jt1->second - data.moment << " "
                << std::sqrt(jt2->second * jt2->second + data.d_moment * data.d_moment) << " ";
        }
        output_file << std::endl;
    }
}
