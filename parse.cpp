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
    // the fly here, and perform averaging later.
    std::map<double, std::vector<double>> energies;
    std::map<double, std::vector<double>> d_energies;
    std::map<double, std::vector<double>> moments;
    std::map<double, std::vector<double>> d_moments;
    std::map<double, std::map<double, std::vector<double>>> Cms;
    std::map<double, std::map<double, std::vector<double>>> d_Cms;
};

struct averaged_data
{
    // It's 2:47 AM and I don't know what I'm writing anymore; there is probably
    // a better way to structure this I will think of later
    std::map<double, double> energy;
    std::map<double, double> d_energy;
    std::map<double, double> moment;
    std::map<double, double> d_moment;
    std::map<double, std::map<double, double>> Cm;
    std::map<double, std::map<double, double>> d_Cm;
};

void parse_file(std::string in_fname, unaveraged_data &data);
averaged_data average_data(unaveraged_data &data);
void print_to_file(std::string out_fname, averaged_data &data);

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
        unaveraged_data data;
        for (auto &p : fs::directory_iterator(fs::current_path()))
        {
            // For now, only examine files from the first realization (*r1.out),
            // and only files that contain the local moment correlations (r*.out)
            // TODO: add support for multiple realizations, and average over them.
            auto str = p.path().filename().string();
            if (str.at(0) == 'r')
            {
                parse_file(str, data);
            }
        }
        auto final_temp_data = average_data(data);
        fs::current_path(fs::current_path().parent_path().parent_path());
        print_to_file("beta_" + inverse_temp + "_parsed.dat", final_temp_data);
    }
    return 0;
}

void parse_file(std::string in_fname, unaveraged_data &data)
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
    data.energies[mu].push_back(avg_energy);
    data.d_energies[mu].push_back(energy_uc);

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
    std::map<double, std::vector<double>> moment_correls, moment_correl_errs;
    std::map<double, double> fin_moment_correls, fin_moment_correl_errs;

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
        data.Cms[mu][r].push_back(cf);
        data.d_Cms[mu][r].push_back(cf_uc);

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

    data.moments[mu].push_back(local_moment);
    data.d_moments[mu].push_back(local_moment_uc);
}

averaged_data average_data(unaveraged_data &data)
{
    averaged_data final_data;

    // Okay, so for a given temperature and chemical potential, we have basically
    // turned a sequence of files with one entry per variable into a single file
    // with a sequence of data entries per variable. Now we have to average over
    // those files (in other words, over these new data entries).

    // So first, iterate over these chemical potentials, since the chemical potential
    // files are stored in the same directory (if they were in separate directories,
    // our lives would be easier)

    // Energy average
    for (auto it1 = data.energies.begin(), it2 = data.d_energies.begin();
            it1 != data.energies.end();
            it1++, it2++)
    {
        const auto mu = it1->first;
        final_data.energy[mu] = 0.;
        final_data.d_energy[mu] = 0.;
        for (auto jt1 = it1->second.begin(), jt2 = it2->second.begin();
                jt1 != it1->second.end();
                jt1++, jt2++)
        {
            final_data.energy[mu] += *jt1;
            final_data.d_energy[mu] += (*jt2) * (*jt2);
        }
        final_data.energy[mu] /= it1->second.size();
        final_data.d_energy[mu] = std::sqrt(final_data.d_energy[mu]) / it1->second.size();
    }

    // Moment average
    for (auto it1 = data.moments.begin(), it2 = data.d_moments.begin();
            it1 != data.moments.end();
            it1++, it2++)
    {
        const auto mu = it1->first;
        final_data.moment[mu] = 0.;
        final_data.d_moment[mu] = 0.;
        for (auto jt1 = it1->second.begin(), jt2 = it2->second.begin();
                jt1 != it1->second.end();
                jt1++, jt2++)
        {
            final_data.moment[mu] += *jt1;
            final_data.d_moment[mu] += (*jt2) * (*jt2);
        }
        final_data.moment[mu] /= it1->second.size();
        final_data.d_moment[mu] = std::sqrt(final_data.d_moment[mu]) / it1->second.size();

        // We care about the squared local moment, so we update the uncertainty
        // (relative errors added in quadrature) and measured values

        // For the uncertainty, recall the following:
        // delta(a*b) = |a*b| sqrt( (delta(a)/a)^2 + (delta(b)/b)^2 )
        // In the case that a=b, this simplifies to
        // delta(a^2) = |a^2| sqrt( 2 (delta(a)/a)^2 )
        //            = |a^2| |delta(a)/a| sqrt(2)
        //            = |a * delta(a)| sqrt(2)
        final_data.d_moment[mu] = std::abs(final_data.moment[mu] * final_data.d_moment[mu]) * std::sqrt(2.);
        final_data.moment[mu] *= final_data.moment[mu];
    }

    // Now, for each radius, we perform averaging of the corr. fun, and
    // add the errors in quadrature

    for (auto it1 = data.Cms.begin(), it2 = data.d_Cms.begin();
            it1 != data.Cms.end();
            it1++, it2++)
    {
        const auto mu = it1->first;
        for (auto jt1 = it1->second.begin(), jt2 = it2->second.begin();
                jt1 != it1->second.end();
                jt1++, jt2++)
        {
            const auto r = jt1->first;
            final_data.Cm[mu][r] = 0;
            final_data.d_Cm[mu][r] = 0;
            for (auto kt1 = jt1->second.begin(), kt2 = jt2->second.begin();
                    kt1 != jt1->second.end();
                    kt1++, kt2++)
            {
                final_data.Cm[mu][r] += *kt1;
                final_data.d_Cm[mu][r] += (*kt2) * (*kt2);
            }
            final_data.Cm[mu][r] /= jt1->second.size();
            final_data.d_Cm[mu][r] = std::sqrt(final_data.d_Cm[mu][r]) / jt1->second.size();
        }
    }

    return final_data;
}

void print_to_file(std::string out_fname, averaged_data &data)
{
    std::ofstream output_file(out_fname);
    // Now, print to files
    output_file << "# mu   <H>    delta_<H>    <mz^2>    delta_<mz^2>    ";
    // Take the first chemical potential, and print out all the lattice distances
    // it doesn't matter which potential we use, since they're all on the same lattice.
    for (auto it = data.Cm.begin()->second.begin(); it != data.Cm.begin()->second.end(); it++)
    {
        auto r = it->first;
        output_file << "Cm(" << std::setprecision(3) << r <<
            ")    dCm(" << std::setprecision(3) << r << ")    ";
    }
    output_file << std::endl;

    // Now, loop over chemical potentials and print the data for each
    for (auto it = data.Cm.begin(); it != data.Cm.end(); it++)
    {
        auto mu = it->first;
        output_file << mu << " " << data.energy[mu] << " " << data.d_energy[mu] << " "
            << data.moment[mu] << " " << data.d_moment[mu] << " ";

        for (auto jt1 = data.Cm[mu].begin(), jt2 = data.d_Cm[mu].begin();
                jt1 != data.Cm[mu].end();
                jt1++, jt2++)
        {
            output_file << std::setprecision(8) << jt1->second - data.moment[mu] << " "
                << std::sqrt(jt2->second * jt2->second + data.d_moment[mu] * data.d_moment[mu]) << " ";
        }
        output_file << std::endl;
    }
}
