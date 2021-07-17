
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "gsl/gsl_spline.h"

int main(int argc, char *argv[])
{
    constexpr auto ninterp = 3u;
    const double beta = std::stod(argv[1]);
    constexpr double dbeta = 0.05;

    constexpr auto nmu = 61u;
    constexpr double dmu = 0.5;

    double n[ninterp][nmu];
    double dn_dt[nmu];

    for (auto i = 0; i < ninterp; i++)
    {
        int lens[ninterp] = {4, 3, 4};
        auto newbeta = beta + (i - 1) * dbeta;
        const auto ifname = std::string("beta_") + std::to_string(newbeta).erase(lens[i], std::string::npos) + "_core.dat";
        std::ifstream in_file(ifname);

        if (!in_file.is_open())
            throw std::runtime_error("Could not open beta file with filename " + ifname);

        std::string line;
        std::getline(in_file, line);

        while (std::getline(in_file, line))
        {
            std::stringstream sstr(line);
            double mu_index;
            sstr >> mu_index;
            mu_index = mu_index * 2 + 30;
            double discard;
            sstr >> discard >> discard >> discard >> discard;
            sstr >> n[i][static_cast<unsigned>(mu_index)];
        }
    }

    std::ofstream dn_dt_file(std::string("dn_dt_beta_") + argv[1] + ".dat");
    std::ofstream dn_dt_fd_file(std::string("dn_dt_beta_") + argv[1] + "_fd.dat");
    // Now, interpolate for each mu
    for (auto i = 0; i < nmu; i++)
    {
        const double xarr[ninterp] = { 1./(beta + dbeta), 1./beta, 1./(beta - dbeta) };
        double yarr[ninterp];

        for (auto j = 0; j < ninterp; j++)
            yarr[j] = n[ninterp - 1 - j][i];

        auto acc = gsl_interp_accel_alloc();
        auto spline = gsl_spline_alloc(gsl_interp_polynomial, ninterp);
        gsl_spline_init(spline, xarr, yarr, ninterp);

        dn_dt[i] = gsl_spline_eval_deriv(spline, 1./beta, acc);
        dn_dt_file << i * 0.5 - 15 << " " << dn_dt[i] << std::endl;
        // std::cerr << dn_dt[i] << std::endl;

        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);

        // Let's try it the other way...
        auto h = 0.5 * (1/(beta - dbeta) - 1/(beta + dbeta));
        // dn_dt[i] = 0.5 * (n[0][i] - n[2][i]) / h;
        dn_dt_fd_file << i * 0.5 - 15 << " " << 0.5 * (n[0][i] - n[2][i]) / h << std::endl;
    }
    dn_dt_file.close();

    // Now, for the integration step. We spline our curve and then integrate, and print out points.

    double xarr[nmu];

    for (auto i = 0; i < nmu; i++)
        xarr[i] = i*0.5 - 15;

    auto acc = gsl_interp_accel_alloc();
    auto spline = gsl_spline_alloc(gsl_interp_cspline, nmu);
    gsl_spline_init(spline, xarr, dn_dt, nmu);

    std::ofstream out_file(std::string("entropy_beta_") + argv[1] + ".dat");
    for (auto i = 0; i < 601; i++)
    {
        const auto x = i * 0.05 - 15;
        const auto s = gsl_spline_eval_integ(spline, -15., x, acc);
        out_file << x << " " << s << std::endl;
    }
    gsl_interp_accel_free(acc);
    gsl_spline_free(spline);

    return 0;
}
