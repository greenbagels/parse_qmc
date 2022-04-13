/*! @file sign_data.cpp
 *  @brief Simple Total Average Sign parsing tool
 *  @author Sameed Pervaiz (pervaiz.8@osu.edu)
 *  @date 2021-09-14
 *  @copyright GPLv3
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>
#include <sstream>

#include "common.hpp"

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "Error: must provide data input directory!\n";
        return 1;
    }

    std::vector<qmc_data> pts;

    namespace fs = std::filesystem;
    for (auto &p : fs::directory_iterator(argv[1]))
    {
        std::string fname = p.path().filename().string();
        if (fname.find("rz") == std::string::npos)
        {
            std::cerr << "Warning: skipping non-data file " << fname << std::endl;
            continue;
        }
        qmc_data data;
        try
        {
            data = parse_qmc_fname(fname);
        }
        catch (std::exception &e)
        {
            std::cerr << "Exception caught in filename parsing: " <<  e.what() << std::endl;
            return 1;
        }
        parse_qmc_file(p.path(), data);
        pts.push_back(data);
    }

    std::cerr << "Parsed " << pts.size() << " data entries\n";

    std::sort(pts.begin(), pts.end(), [&](qmc_data &a, qmc_data &b) {return (a.mu < b.mu) || ((a.mu == b.mu) && (a.beta < b.beta));});

    std::ofstream output_file("total_signs.dat");
    output_file << "#mu    beta    <total sign>\n";
    for (auto i = pts.begin(); i != pts.end(); i++)
    {
        output_file << i->mu << " " << i->density << " " << i->beta << " " << i->average_total_sign << std::endl;
        if (i+1 != pts.end())
        {
            if (i->mu != (i+1)->mu)
                output_file << std::endl;
        }
    }


    return 0;
}

