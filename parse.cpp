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
        std::ofstream parsed_file("beta_" + inverse_temp + "_parsed.dat");
        fs::current_path(p);
        for (auto &p : fs::directory_iterator(fs::current_path()))
        {
            // For now, only examine files from the first realization (*r1.out),
            // and only files that contain the local moment correlations (r*.out)
            auto str = p.path().filename().string();
            if (str.at(0) == 'r' && str.substr(str.size()-6) == "r1.out")
            {
                std::cerr << str << std::endl;
                // Now, extract the chemical potential, I guess
                std::stringstream sstream(str.substr(str.find("mu")+2));
                double mu;
                sstream >> mu;
                std::cout << "Parsing chemical potential \u03BC = " << mu << std::endl;
                std::ifstream data_file(str);
                std::ofstream correl_r_file("cf_beta_" + inverse_temp + "_mu_" + std::to_string(mu) + ".dat");

                // Now, parse the file. For now, let's only track the average site energy
                // So skip everything up until that line
                std::string line;
                do
                {
                    std::getline(data_file, line);
                }
                while (line.find("Average Energy") == line.npos);
                std::cout << line << std::endl;
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

                    // First, construct the vector if empty
                    if (moment_correls.find(r) == moment_correls.end())
                    {
                        // TODO: there *has* to be a better way to do this
                        moment_correls[r] = std::vector<double>();
                        moment_correl_errs[r] = std::vector<double>();
                    }
                    // Now, tack on our values at the given radius
                    moment_correls.at(r).push_back(cf);
                    moment_correl_errs.at(r).push_back(cf_uc);

                    std::getline(data_file, line);
                }

                // Now, for each radius, we perform averaging of the corr. fun, and
                // add the errors in quadrature

                for (auto it1 = moment_correls.begin(), it2 = moment_correl_errs.begin();
                        it1 != moment_correls.end();
                        it1++, it2++)
                {
                    const auto r = it1->first;
                    fin_moment_correls[r] = 0.;
                    fin_moment_correl_errs[r] = 0.;
                    for (auto i1 = it1->second.begin(), i2 = it2->second.begin();
                            i1 != it1->second.end();
                            i1++, i2++)
                    {
                        fin_moment_correls[r] += *i1;
                        fin_moment_correl_errs[r] += (*i2) * (*i2);
                    }
                    fin_moment_correls[r] /= it1->second.size();
                    fin_moment_correl_errs[r] = std::sqrt(fin_moment_correl_errs[r]);
                }

                /*
                // Now, discard the (0,0) line, AKA the variance
                std::getline(data_file, line);

                // Now, the nearest-neighbor (0,1) line
                std::getline(data_file, line);
                sstream = std::stringstream(line);

                double n_moment, n_moment_uc;
                // X Y MOMENT MOMENT_UNCERTAINTY
                sstream >> n_moment >> n_moment >> n_moment >> n_moment_uc;

                // Now, the next-nearest-neighbor (1,1) line
                do
                {
                    std::getline(data_file, line);
                }
                while (line.find("1   1") == line.npos);
                sstream = std::stringstream(line);
                double nn_moment, nn_moment_uc;
                sstream >> nn_moment >> nn_moment >> nn_moment >> nn_moment_uc;
                */

                // Now, local moment.
                do
                {
                    std::getline(data_file, line);
                }
                while (line.find("0 0 local moment zz=") == line.npos);

                sstream = std::stringstream(line.substr(line.find("=")+1));
                double local_moment, local_moment_uc;
                sstream >> local_moment >> local_moment_uc;

                // We care about the squared local moment, so we update the uncertainty
                // (relative errors added in quadrature) and measured values

                // For the uncertainty, recall the following:
                // delta(a*b) = |a*b| sqrt( (delta(a)/a)^2 + (delta(b)/b)^2 )
                // In the case that a=b, this simplifies to
                // delta(a^2) = |a^2| sqrt( 2 (delta(a)/a)^2 )
                //            = |a^2| |delta(a)/a| sqrt(2)
                //            = |a * delta(a)| sqrt(2)
                local_moment_uc = std::abs(local_moment * local_moment_uc) * std::sqrt(2.);
                local_moment *= local_moment;

                // Now, print to files
                parsed_file << "# mu   <H>    delta_<H>    <mz^2>    delta_<mz^2>    ";
                for (auto i = fin_moment_correls.begin(); i != fin_moment_correls.end(); i++)
                {
                    auto r = i->first;
                    parsed_file << "Cm(" << std::setprecision(3) << r <<
                        ")    dCm(" << std::setprecision(3) << r << ")    ";
                }
                parsed_file << std::endl;

                parsed_file << mu << " " << avg_energy << " " << energy_uc << " "
                    << local_moment << " " << local_moment_uc << " ";

                for (auto i1 = fin_moment_correls.begin(), i2 = fin_moment_correl_errs.begin();
                        i1 != fin_moment_correls.end();
                        i1++, i2++)
                {
                    parsed_file << std::setprecision(8) << i1->second - local_moment << " "
                        << std::sqrt(i2->second * i2->second + local_moment_uc * local_moment_uc) << " ";
                }
                parsed_file << std::endl;
            }
        }
        fs::current_path(fs::current_path().parent_path().parent_path());
    }
    return 0;
}
