
using KahanSummation
using Measurements

include("sample.jl")

function parse_bandmott_header(filename)
    m = match(r"rz(\d+)l(\d+)u(\d+\.\d+)dt(0\.\d+)mu(\d+\.\d+)(n?)r(\d+)\.out", filename)
    if m.captures[6] == "n"
        return Header(parse(Int, m[1]),
                      parse(Int, m[2]),
                      parse(Float64, m[3]),
                      parse(Float64, m[4]),
                     -parse(Float64, m[5]),
                      parse(Int, m[2]) * parse(Float64, m[4]))
    end

    return Header(parse(Int, m[1]),
                  parse(Int, m[2]),
                  parse(Float64, m[3]),
                  parse(Float64, m[4]),
                  parse(Float64, m[5]),
                  parse(Int, m[2]) * parse(Float64, m[4]))
end

function parse_bandmott_line(target_str, lines)
    i = findfirst(x -> occursin(target_str, x), lines)
    line = lines[i]
    # Let's hope that popping the front makes future parses faster
    for j in 1:i
        popfirst!(lines)
    end
    strvals = line[findfirst(isequal('='), line)+1:end]
    vals = [parse(Float64, str) for str in split(strvals)]
    if vals[2] == 0.
        return [vals[1], 1e-8]
    end
    return vals
end

function parse_bandmott_matrix(target_str, lines, size)
    n = Int(size/2 + 1)
    mat = Matrix{Measurement{Float64}}(undef, n, n)

    i = findfirst(x -> occursin(target_str, x), lines) + 1
    for line in lines[i:end]
        str = split(line)
        x = tryparse(Int, str[1])
        y = tryparse(Int, str[2])
        if isnothing(x) || isnothing(y)
            return mat
        end
        if x > y || isnothing(tryparse(Float64, str[3]))
            continue
        end
        x += 1
        y += 1
        val = 0
        uc = 0
        if str[4] == "+-"
            # This really only happens with the density correlator, for which
            # the situation kinda looks like this:
            # x y nup-nup +- error  nup-ndn +- error

            # By symmetry, nup-nup = ndn-ndn. And so the full density correlator
            # is just 2 up-(up+dn) = 2(up-up + up-dn)
            mat[y, x] = 2 * measurement(parse(Float64, str[3]) + parse(Float64, str[6]),
                                        parse(Float64, str[5]) + parse(Float64, str[8]))
        else
            # This is just a normal line that looks like
            # x y val uncertainty
            mat[y, x] = measurement(parse(Float64, str[3]), parse(Float64, str[4]))
        end
    end
end

function parse_bandmott_file(filename)
    file = open(filename)
    lines = readlines(file)

    vals    = parse_bandmott_line("Average total sign", lines)
    total_sign  = measurement(vals[1], vals[2])

    vals    = parse_bandmott_line("Average density", lines)
    density = measurement(vals[1], vals[2])

    vals    = parse_bandmott_line("Average up occupancy", lines)
    nup     = measurement(vals[1], vals[2])

    vals    = parse_bandmott_line("Average dn occupancy", lines)
    ndn     = measurement(vals[1], vals[2])

    vals    = parse_bandmott_line("Average Energy", lines)
    energy  = measurement(vals[1], vals[2])

    vals    = parse_bandmott_line("Average Nup*Ndn", lines)
    nupndn  = measurement(vals[1], vals[2])

    vals    = parse_bandmott_line("AF correlation function (zz)", lines)
    afsf    = measurement(vals[1], vals[2])

    vals    = parse_bandmott_line("Ferro corr. func. (zz)", lines)
    fsf     = measurement(vals[1], vals[2])

    ninj = parse_bandmott_matrix("density-density correlation fn", lines, 10)
    didj = parse_bandmott_matrix("nud-nud correlation function", lines, 10)
    mimj = parse_bandmott_matrix("mi2x-mi2x correlation function", lines, 10)
    nidj = parse_bandmott_matrix("nup-nud correlation function", lines, 10)
    sisj = parse_bandmott_matrix("zz Spin correlation function", lines, 10)

    close(file)
    return Sample(total_sign, energy, density, nup, ndn, nupndn, measurement(0),
                  measurement(0), fsf, afsf, [ninj, nidj, didj, mimj, sisj])
end

function average_data(unaveraged_data)
    sample_map = Dict{Float64, Sample{Float64}}()
    for (mu, arr) in pairs(unaveraged_data)
        total_sign = weightedmean([sample.total_sign for sample in arr])
        avg_energy = weightedmean([sample.avg_energy for sample in arr])
        n = weightedmean([sample.n for sample in arr])
        nup = weightedmean([sample.nup for sample in arr])
        ndn = weightedmean([sample.ndn for sample in arr])
        nupndn = weightedmean([sample.nupndn for sample in arr])
        z_spin = (nup - ndn) / 2
        moment = n - 2 * nupndn
        ferro_cf = weightedmean([sample.ferro_cf for sample in arr])
        antiferro_cf = weightedmean([sample.antiferro_cf for sample in arr])

        connected = [n^2, n * nupndn, nupndn^2, moment^2, z_spin^2]

        cfs = Array{Matrix{Measurement{Float64}}, 1}(undef, 5)
        for i in 1:5
            width = 2*(size(arr[1].cfs[1], 1)-1)
            cfs[i] = Matrix{Measurement{Float64}}(undef, width, width)
            for j in 1:size(arr[1].cfs[1], 1), k in j:size(arr[1].cfs[1], 1)
                cfs[i][k, j] = weightedmean([sample.cfs[i][k, j] for sample in arr]) - connected[i]
                cfs[i][j, k] = cfs[i][k, j]
            end
        end
        sample_map[mu] = Sample(total_sign, avg_energy, n, nup, ndn, nupndn, z_spin, moment,
                                ferro_cf, antiferro_cf, cfs)
    end
    return sample_map
end

function flesh_out_cfs!(averaged_data)
    for (mu, sample) in averaged_data
        n = Int(size(sample.cfs[1], 1)/2)
        fulln = size(sample.cfs[1], 1)
        for i in 1:5, x in 0:n, y in 0:n
            sample.cfs[i][mod(-y, fulln) + 1, x + 1] =
                sample.cfs[i][y + 1, x + 1]
            sample.cfs[i][mod(-y, fulln) + 1, mod(-x, fulln) + 1] =
                sample.cfs[i][y + 1, x + 1]
            sample.cfs[i][y + 1, mod(-x, fulln) + 1] =
                sample.cfs[i][y + 1, x + 1]
        end
    end
end

function calculate_SAFs(averaged_data)
    # (pi,pi) fourier-transformed spin correlation function
    SAFs = Dict{Float64, Measurement{Float64}}()
    qx = pi
    qy = pi
    for (mu, sample) in averaged_data
        SAFs[mu] = 0
        spin_cf = sample.cfs[5]
        #display(spin_cf)
        #println()
        for ix in 1:size(spin_cf, 1), iy in 1:size(spin_cf, 1),
            jx in 1:size(spin_cf, 1), jy in 1:size(spin_cf, 1)
            dx = ix - jx
            dy = iy - jy
            absdx = Int(abs(dx))
            absdy = Int(abs(dy))
            SAFs[mu] += cos(qx * dx + qy * dy) * spin_cf[absdy + 1, absdx + 1]
        end
        SAFs[mu] /= 100
    end
    return SAFs
end

