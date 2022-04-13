
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>

#include "gsl/gsl_spline.h"

template <typename T>
struct mu_sample
{
    mu_sample(T mu_ = 0., T val_ = 0.) : mu(mu_), val(val_) {}
    T mu;
    T val;
};

template <size_t N>
struct spline_data
{
    double mu;
    std::array<double, N> xarr;
    std::array<double, N> yarr;
    gsl_interp_accel *accelerator;
    gsl_spline *spline;
};

template <typename T>
std::array<T,2> simpsons_rule(T *f, T *df, T h, unsigned npts);

int main(int argc, char *argv[])
{
    constexpr auto ninterp = 5u;
    const double beta = std::stod(argv[1]);
    constexpr double dbeta = 0.05;

    // ninterp is ordered from low beta to high beta
    // so in other words, inverse-temperature ordering!
    std::array<std::vector<mu_sample<double>>, ninterp> n;
    std::array<std::vector<mu_sample<double>>, ninterp> n_err;
    std::array<std::vector<mu_sample<double>>, ninterp> mz;
    std::vector<mu_sample<double>> dn_dt;
    std::vector<mu_sample<double>> symmetrized_dn_dt;
    /*
    std::vector<mu_sample<double>> dn_dt_fd;
    std::vector<mu_sample<double>> symmetrized_dn_dt_fd;
    std::vector<mu_sample<double>> dn_dt_fd_cr;
    std::vector<mu_sample<double>> symmetrized_dn_dt_fd_cr;
    */
    std::vector<mu_sample<double>> dn_dt_err;
    std::vector<mu_sample<double>> symmetrized_dn_dt_err;

    for (auto i = 0; i < ninterp; i++)
    {
        int lens[ninterp] = {4, 4, 3, 4, 4};
        auto newbeta = beta + (i - 2) * dbeta;
        const auto ifname = std::string("beta_") + std::to_string(newbeta).erase(lens[i], std::string::npos) + "_core.dat";
        std::ifstream in_file(ifname);

        if (!in_file.is_open())
            throw std::runtime_error("Could not open beta file with filename " + ifname);

        std::string line;
        std::getline(in_file, line);

        while (std::getline(in_file, line))
        {
            std::stringstream sstr(line);
            mu_sample<double> sample;
            sstr >> sample.mu;
            for (auto j = 0; j < 3; j++)
                sstr >> sample.val; // The last iteration contains the actual value of n
            mz[i].push_back(sample);
            for (auto j = 0; j < 2; j++)
                sstr >> sample.val; // The last iteration contains the actual value of n
            n[i].push_back(sample);
            sstr >> sample.val;
            n_err[i].push_back(sample);
        }
        // Sort on mu, since grid spacing could be non-constant.
        std::sort(mz[i].begin(), mz[i].end(), [](mu_sample<double> a, mu_sample<double> b) { return a.mu < b.mu; } );
        std::sort(n[i].begin(), n[i].end(), [](mu_sample<double> a, mu_sample<double> b) { return a.mu < b.mu; } );
        std::sort(n_err[i].begin(), n_err[i].end(), [](mu_sample<double> a, mu_sample<double> b) { return a.mu < b.mu; } );
    }

    // Run-time check: if one of the n vectors has different mu values from another, we have a problem.
    for (auto it = n.begin(); it != n.end() - 1; it++)
    {
        auto vec1 = *it;
        auto vec2 = *(it+1);
        auto diff = static_cast<int>(it - n.begin());

        // First, check if the samples have the same size, to avoid range errors
        if (vec1.size() != vec2.size())
            throw std::runtime_error("Input files (i=" + std::to_string(diff) + ") do not contain the same grid spacing!");

        // Now, check the values actually match
        for (auto i = 0; i < vec1.size(); i++)
        {
            if (vec1[i].mu != vec2[i].mu)
                throw std::runtime_error("Input files (i=" + std::to_string(diff) + ") do not contain the same grid spacing!");
        }
    }

    std::ofstream dn_dt_file(std::string("dn_dt_beta_") + argv[1] + ".dat");
    // std::ofstream dn_dt_fd_file(std::string("dn_dt_beta_") + argv[1] + "_fd.dat");
    // std::ofstream dn_dt_fd_cr_file(std::string("dn_dt_beta_") + argv[1] + "_fd_cr.dat");
    std::ofstream sym_dn_dt_file(std::string("sym_dn_dt_beta_") + argv[1] + ".dat");
    // std::ofstream sym_dn_dt_fd_file(std::string("sym_dn_dt_beta_") + argv[1] + "_fd.dat");
    // std::ofstream sym_dn_dt_fd_cr_file(std::string("sym_dn_dt_beta_") + argv[1] + "_fd_cr.dat");

    // First, interpolate along each isochime and store the derivatives
    std::vector<spline_data<ninterp>> three_pt_splines;
    // TODO: check for empty vecs
    dn_dt_file << "# mu    mz    n    dn/dt    abs.err\n";
    for (auto i = 0; i < n[0].size(); i++)
    {
        /*
        spline_data<ninterp> lagrange;
        lagrange.mu = n[0][i].mu;
        lagrange.xarr = { 1./(beta + dbeta), 1./beta, 1./(beta - dbeta) };

        for (auto j = 0; j < ninterp; j++)
            lagrange.yarr[j] = n[ninterp - 1 - j][i].val;

        lagrange.accelerator = gsl_interp_accel_alloc();
        lagrange.spline = gsl_spline_alloc(gsl_interp_polynomial, ninterp);
        gsl_spline_init(lagrange.spline, lagrange.xarr.data(), lagrange.yarr.data(), ninterp);
        three_pt_splines.push_back(lagrange);

        mu_sample<double> sample;
        sample.mu = n[0][i].mu;
        sample.val = gsl_spline_eval_deriv(lagrange.spline, 1./beta, lagrange.accelerator);
        dn_dt.push_back(sample);
        dn_dt_file << sample.mu << " " << sample.val << std::endl;
        // std::cerr << dn_dt[i] << std::endl;

        // gsl_spline_free(spline);
        // gsl_interp_accel_free(acc);

        // Let's try it the other way... (naive almost-finite-difference)
        auto h = 0.5 * (1/(beta - dbeta) - 1/(beta + dbeta));
        // dn_dt[i] = 0.5 * (n[0][i] - n[2][i]) / h;
        // TODO: avoid magic numbers 0,2
        sample.val = 0.5 * (n[0][i].val - n[2][i].val) / h;
        dn_dt_fd.push_back(sample);
        dn_dt_fd_file << sample.mu << " " << sample.val << std::endl;
        */
        // Now the chain-rule way, because we were being dumb earlier!
        auto h = dbeta;
        mu_sample<double> sample;
        sample.mu = n[0][i].mu;
        sample.val = (8.*(n[3][i].val - n[1][i].val) + n[0][i].val - n[4][i].val) / (12. * h);
        // Tack on the -1/T^2 = - beta^2
        sample.val *= - beta * beta;

        dn_dt.push_back(sample);
        dn_dt_file << sample.mu << " " << mz[2][i].val << " " << n[2][i].val << " " << sample.val;

        // So now we account for the uncertainty; take {beta, h} to be exact,
        // and propagate the errors from the two n values; this essentially means
        // the error is beta^2 / (2h) * sqrt(dn1^2 + dn2^2)
        sample.val = beta * beta / (12. * h)
            * std::hypot(8. * std::hypot(n_err[3][i].val, n_err[1][i].val), std::hypot(n_err[0][i].val, n_err[4][i].val));
        dn_dt_err.push_back(sample);
        dn_dt_file << " " << sample.val << std::endl;
    }

    /*
    std::ofstream surface_file(std::string("density_surface_beta_") + argv[1] + ".dat");
    // Now, interpolate along mu for each T value
    const auto n_t_pts = 50;
    for (auto i = 0; i <= n_t_pts; i++)
    {
        std::vector<double> xarr, yarr;

        // First, pick the T value that we interpolate along (with respect to mu)
        auto T_start = 1./(beta+dbeta);
        auto T_end = 1./(beta-dbeta);
        auto dT = (T_end - T_start) / n_t_pts;
        auto T = T_start + i * dT;

        // Now, build the list of n points as a function of mu.
        for (auto j = three_pt_splines.begin(); j != three_pt_splines.end(); j++)
        {
            xarr.push_back(j->mu);
            yarr.push_back(gsl_spline_eval(j->spline, T, j->accelerator));
        }

        // Now, create the spline! Use cubic splines to avoid bad endpoint behavior...
        auto accel = gsl_interp_accel_alloc();
        auto spline = gsl_spline_alloc(gsl_interp_cspline, xarr.size());
        gsl_spline_init(spline, xarr.data(), yarr.data(), xarr.size());

        auto n_mu_pts = 150;
        for (auto j = 0; j <= n_mu_pts; j++)
        {
            auto mu = -15. + j*(30./n_mu_pts);
            auto n_surf_pt = gsl_spline_eval(spline, mu, accel);
            surface_file << T << " " << mu << " " << n_surf_pt << std::endl;
        }
        surface_file << std::endl;

        gsl_interp_accel_free(accel);
        gsl_spline_free(spline);

    }
    */
    // dn_dt_file.close();

    // Now, for the integration step. We spline our curve and then integrate, and print out points.

    // But first, symmetrize the data

    sym_dn_dt_file << "# mu    mz    n    dn/dt    abs.err\n";
    for (auto i = 0; i < dn_dt.size(); i++)
    {
        auto reverse_i = dn_dt.size() - (i + 1);
        mu_sample<double> sample;
        sample.mu = dn_dt[i].mu;

        /*
        sample.mu = dn_dt[i].mu;
        sample.val = 0.5 * (dn_dt[i].val - dn_dt[reverse_i].val);
        sym_dn_dt_file << sample.mu << " " << sample.val << std::endl;
        symmetrized_dn_dt.push_back(sample);

        // Re-use the same mu
        sample.val = 0.5 * (dn_dt_fd[i].val - dn_dt_fd[reverse_i].val);
        sym_dn_dt_fd_file << sample.mu << " " << sample.val << std::endl;
        symmetrized_dn_dt_fd.push_back(sample);
        */

        // Re-use the same mu
        sample.val = 0.5 * (dn_dt[i].val - dn_dt[reverse_i].val);
        sym_dn_dt_file << sample.mu << " " << mz[2][i].val << " " << n[2][i].val << " " << sample.val;
        symmetrized_dn_dt.push_back(sample);

        // Now, the error; treat 0.5 as exact, add the other erros in quadrature
        sample.val = 0.5 * std::hypot(dn_dt_err[i].val, dn_dt_err[reverse_i].val);
        sym_dn_dt_file << " " << sample.val << std::endl;
        symmetrized_dn_dt_err.push_back(sample);
    }

    /*
    std::vector<double> xarr;
    std::vector<double> yarr;

    for (auto i = 0; i < symmetrized_dn_dt.size(); i++)
    {
        xarr.push_back(symmetrized_dn_dt[i].mu);
        yarr.push_back(symmetrized_dn_dt[i].val);
    }

    auto acc = gsl_interp_accel_alloc();
    auto spline = gsl_spline_alloc(gsl_interp_cspline, xarr.size());
    gsl_spline_init(spline, xarr.data(), yarr.data(), xarr.size());

    std::vector<double> muarr, narr;

    for (auto i = 0; i < n[2].size(); i++)
    {
        muarr.push_back(n[2][i].mu);
        narr.push_back(n[2][i].val);
    }

    auto n_acc = gsl_interp_accel_alloc();
    auto n_spline = gsl_spline_alloc(gsl_interp_cspline, muarr.size());
    gsl_spline_init(n_spline, muarr.data(), narr.data(), muarr.size());

    std::ofstream out_file(std::string("entropy_beta_") + argv[1] + ".dat");
    for (auto i = 0; i <= 600; i++)
    {
        const auto mu = i * 0.05 - 15.;
        const auto density = gsl_spline_eval(n_spline, mu, n_acc);
        const auto s = gsl_spline_eval_integ(spline, -15., mu, acc);
        out_file << density << " " << s << std::endl;
    }
    gsl_interp_accel_free(acc);
    gsl_spline_free(spline);
    gsl_interp_accel_free(n_acc);
    gsl_spline_free(n_spline);
    */

    // First, integrate the left segment, [xmin = -15, xmax = -10, dx = 0.5]
    // Which means n = 11

    // If n=1, then we use the trapezoidal rule. If n>1, then we split into the
    // two cases of even n, or odd n.

    std::vector<mu_sample<double>> integral;
    std::vector<mu_sample<double>> integral_err;
    std::vector<double> samples, dsamples;
    /*
    for (auto i = 0; i <= 10; i++)
    {
        samples.push_back(symmetrized_dn_dt[i].val);
        dsamples.push_back(symmetrized_dn_dt_err[i].val);
        auto result = simpsons_rule(samples.data(), dsamples.data(), 0.5, i+1);
        integral.emplace_back(-15. + i * 0.5, result[0]);
        integral_err.emplace_back(-15. + i * 0.5, result[1]);
    }
    for (auto i = 1; i <= 100; i++)
    {
        auto offset = 10;
        samples.push_back(symmetrized_dn_dt[i+offset].val);
        dsamples.push_back(symmetrized_dn_dt_err[i+offset].val);
        auto result = simpsons_rule(samples.data()+offset, dsamples.data()+offset, 0.1, i+1);
        integral.emplace_back(-10. + i * 0.1, result[0] + integral[offset].val);
        integral_err.emplace_back(-10. + i * 0.1, result[1] + integral_err[offset].val);
    }
    */
    /* // This is for when we ignore forced symmetry
    for (auto i = 1; i <= 10; i++)
    {
        auto offset = 210;
        samples.push_back(symmetrized_dn_dt[i+offset].val);
        dsamples.push_back(symmetrized_dn_dt_err[i+offset].val);
        auto result = simpsons_rule(samples.data()+offset, dsamples.data()+offset, 0.5, i+1);
        integral.emplace_back(10. + i * 0.5, result[0] + integral[offset].val);
        integral_err.emplace_back(10. + i * 0.5, result[1] + integral_err[offset].val);
    }*/
    /*
    for (auto i = 1; i <= 110; i++)
    {
        auto end = 110;
        integral.emplace_back(-integral[end - i].mu, integral[end - i].val);
        integral_err.emplace_back(-integral_err[end - i].mu, integral_err[end - i].val);
    }
    */
    // dmu=0.2 integral
    for (auto i = 0; i <= 75; i++)
    {
        samples.push_back(symmetrized_dn_dt[i].val);
        dsamples.push_back(symmetrized_dn_dt_err[i].val);
        auto result = simpsons_rule(samples.data(), dsamples.data(), 0.2, i+1);
        integral.emplace_back(-15. + i * 0.2, result[0]);
        integral_err.emplace_back(-15. + i * 0.2, result[1]);
    }
    // symmetrize
    for (auto i = 1; i <= 75; i++)
    {
        auto end = 75;
        integral.emplace_back(-integral[end - i].mu, integral[end - i].val);
        integral_err.emplace_back(-integral_err[end - i].mu, integral_err[end - i].val);
    }

    std::ofstream entropy_simpson(std::string("entropy_simpson_beta_") + argv[1] + ".dat");
    entropy_simpson << "# mu    n    mz^2    S    ds\n";

    for (auto i = 0; i < integral.size(); i++)
    {
        const auto mu = integral[i].mu;
        const auto density = n[2][i].val;
        const auto moment = mz[2][i].val;
        entropy_simpson << mu << " " << density << " " << moment << " "
            << integral[i].val << " " << integral_err[i].val << std::endl;
    }
    return 0;
}


// Integrates a function using a specified sampling of its points stored in the
// pointer f
template <typename T>
std::array<T,2> simpsons_rule(T *f, T *df, T h, unsigned npts)
{
    auto nintervals = npts - 1;
    std::array<T,2> integ;

    if (nintervals == 0)
    {
        integ[0] = 0.;
        integ[1] = 0.;
    }
    else if (nintervals == 1)
    {
        // Do trapezoidal rule
        integ[0] = 0.5 * h * (f[0] + f[1]);
        integ[1] = 0.5 * h * std::hypot(df[0], df[1]);
    }
    else if (nintervals % 3 == 0)
    {
        // Do simpson's 3/8 rule
        auto end = npts - 1;
        integ[0] = f[0] + f[end];
        integ[1] = df[0] * df[0] + df[end] * df[end];
        double weights[3] = {2, 3, 3};
        for (auto j = 1; j < end; j++)
        {
            integ[0] += weights[j % 3] * f[j];
            auto errterm = weights[j % 3] * df[j];
            integ[1] += errterm * errterm;
        }
        integ[0] *= 3./8. * h;
        integ[1] = 3./8. * h * std::sqrt(integ[1]);
    }
    else
    {
        // Do as many 3/8 rules as possible, and then do 1/3 rules for the remainder
        int n3, n2;
        switch (nintervals % 3)
        {
            // 3k+1 == 3(k-1)+4, which is of the form 3k+2j
            case 1:
                n3 = nintervals / 3 - 1;
                n2 = 2;
                break;
            // 3k+2 is already of the form 3k+2j
            case 2:
                n3 = nintervals / 3;
                n2 = 1;
                break;
            // This'll never happen, as we handled it above already
            default:
                throw std::runtime_error("nintervals is both divisible and indivisible by 3?\n");
                break;
        }
        auto sub1 = 0., sub2 = 0., sub1_err = 0., sub2_err = 0.;
        auto end = 0;
        if (n3 > 0)
        {
            end = n3 * 3;
            sub1 = f[0] + f[end];
            sub1_err = df[0] * df[0] + df[end] * df[end];
            double weights[3] = {2, 3, 3};
            for (auto j = 1; j < end; j++)
            {
                sub1 += weights[j % 3] * f[j];
                auto errterm = weights[j % 3] * df[j];
                sub1_err += errterm * errterm;
            }
            sub1 *= 3./8. * h;
            sub1_err *= 9./64. * h * h;
        }

        // Now do the 1/3 rule.
        if (n2 > 0)
        {
            f = f + end;
            df = df + end;
            end = n2 * 2;
            sub2 = f[0] + f[end];
            sub2_err = df[0] * df[0] + df[end] * df[end];
            double weights[2] = {2, 4};
            for (auto j = 1; j < end; j++)
            {
                sub2 += weights[j % 2] * f[j];
                auto errterm = weights[j % 2] * df[j];
                sub2_err += errterm * errterm;
            }
            sub2 *= h/3.;
            sub2_err *= h*h/9.;
        }
        integ[0] = sub1 + sub2;
        integ[1] = std::sqrt(sub1_err + sub2_err);
    }
    return integ;
}
