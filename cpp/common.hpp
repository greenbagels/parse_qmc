#ifndef QMC_COMMON_HPP
#define QMC_COMMON_HPP

#include <string>
#include "quantity.hpp"

struct qmc_data
{
    // These are fixed parameters that are known BEFORE the simulation starts
    unsigned n;
    unsigned l;
    unsigned u;
    double mu;
    double beta;
    double dt;
    // This is determined after running the sim
    double average_total_sign;
    double density;
};

template <typename FLOAT>
class matrix
{
    public:
        matrix(std::size_t n)
        {
            data = std::make_unique<qmc::quantity<FLOAT>[]>(n*n);
            this->n = n;
        }
        matrix(matrix&& mat) : n(mat.n), data(std::move(mat.data))
        {
        }
        matrix& operator=(matrix&& mat)
        {
            n = mat.n;
            data = std::move(mat.data);
            return *this;
        }

        qmc::quantity<FLOAT> &operator()(std::size_t x, std::size_t y)
        {
            if (x > this->width() || y > this->width() || x < 0 || y < 0)
            {
                throw std::out_of_range("Matrix subscript(s) out of range!\n");
            }
            return data[y * this->width() + x];
        }

        const std::size_t width()
        {
            return n;
        }

    private:
        std::size_t n;
        std::unique_ptr<qmc::quantity<FLOAT>[]> data;
};

qmc_data parse_qmc_fname(std::string fname)
{
    qmc_data data;
    // The file in general will look like rz[INT]l[INT]u[FLOAT]dt[FLOAT]mu[FLOAT](n)r[INT].out
    // The first thing we should check is whether to expect nr or just r.
    // TODO: error-check, loop this
    auto pos = fname.find('l');
    auto prev_pos = 0;
    data.n = std::stoi(fname.substr(prev_pos+2, pos));
    prev_pos = pos;

    pos = fname.find('u', prev_pos);
    data.l = std::stoi(fname.substr(prev_pos+1, pos));
    prev_pos = pos;

    pos = fname.find("dt", prev_pos);
    data.u = std::stod(fname.substr(prev_pos+1, pos));
    prev_pos = pos;

    pos = fname.find("mu", prev_pos);
    data.dt = std::stod(fname.substr(prev_pos+2, pos));
    prev_pos = pos;

    pos = fname.find("r", prev_pos);
    if (fname[pos-1] == 'n')
    {
        // Don't edit pos; if we ever want to parse r, we can rely on pos+1 being the start
        // of the integer.
        data.mu = -std::stod(fname.substr(prev_pos+2, pos-1));
    }
    else
    {
        data.mu = std::stod(fname.substr(prev_pos+2, pos));
    }

    // If desired, parse r here. (not implemented)

    // Now, just calculate beta.
    data.beta = data.l * data.dt;
    return data;
}

template <typename T>
int parse_qmc_file(T qmc_file_handle, qmc_data &data)
{
    std::ifstream qmc_file(qmc_file_handle);
    if (!qmc_file.is_open())
    {
        std::cerr << "Error: could not open file\n";
        return 1;
    }

    std::string line;
    while (std::getline(qmc_file, line))
    {
        if (line.find("Average total sign") != std::string::npos)
        {
            line = line.substr(line.find("=")+1);
            std::stringstream sstr(line);
            sstr >> data.average_total_sign;

            std::getline(qmc_file, line);
            line = line.substr(line.find("=")+1);
            sstr = std::stringstream(line);
            sstr >> data.density;
            return 0;
        }
    }

    data.average_total_sign = std::nan("");
    return 1;
}

#endif
