
using KahanSummation
using Measurements

struct Header
    n::Int
    l::Int
    u::Float64
    dt::Float64
    mu::Float64
    beta::Float64
end

struct Datum
    total_sign::Measurement{Float64}
    avg_energy::Measurement{Float64}
    n::Measurement{Float64}
    nup::Measurement{Float64}
    ndn::Measurement{Float64}
    nupndn::Measurement{Float64}
    z_spin::Measurement{Float64}
    moment::Measurement{Float64}
    ferro_cf::Measurement{Float64}
    antiferro_cf::Measurement{Float64}
    cfs::Array{Matrix{Measurement{Float64}}, 1}
end

function parse_header(filename)
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

function parse_line(target_str, lines)
    i = findfirst(x -> occursin(target_str, x), lines)
    line = lines[i]
    # Let's hope that popping the front makes future parses faster
    for j in 1:i
        popfirst!(lines)
    end
    strvals = line[findfirst(isequal('='), line)+1:end]
    vals = [parse(Float64, str) for str in split(strvals)]
    if vals[2] == 0.
        return [vals[1], 1e-7]
    end
    return vals
end

function parse_matrix(target_str, lines, size)
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
            mat[y, x] = 2*measurement(parse(Float64, str[3]) + parse(Float64, str[6]),
                                      parse(Float64, str[5]) + parse(Float64, str[8]))
        else
            # This is just a normal line that looks like
            # x y val uncertainty
            mat[y, x] = measurement(parse(Float64, str[3]), parse(Float64, str[4]))
        end
    end
end

function parse_file(filename)
    file = open(filename)
    lines = readlines(file)

    vals    = parse_line("Average total sign", lines)
    total_sign  = measurement(vals[1], vals[2])

    vals    = parse_line("Average density", lines)
    density = measurement(vals[1], vals[2])

    vals    = parse_line("Average up occupancy", lines)
    nup     = measurement(vals[1], vals[2])

    vals    = parse_line("Average dn occupancy", lines)
    ndn     = measurement(vals[1], vals[2])

    vals    = parse_line("Average Energy", lines)
    energy  = measurement(vals[1], vals[2])

    vals    = parse_line("Average Nup*Ndn", lines)
    nupndn  = measurement(vals[1], vals[2])

    vals    = parse_line("AF correlation function (zz)", lines)
    afsf    = measurement(vals[1], vals[2])

    vals    = parse_line("Ferro corr. func. (zz)", lines)
    fsf     = measurement(vals[1], vals[2])

    ninj = parse_matrix("density-density correlation fn", lines, 10)
    didj = parse_matrix("nud-nud correlation function", lines, 10)
    mimj = parse_matrix("mi2x-mi2x correlation function", lines, 10)
    nidj = parse_matrix("nup-nud correlation function", lines, 10)
    sisj = parse_matrix("zz Spin correlation function", lines, 10)

    close(file)
    return Datum(total_sign, energy, density, nup, ndn, nupndn, 0, 0,
                 fsf, afsf, [ninj, nidj, didj, mimj, sisj])
end

function average_data(unaveraged_data)
    datum_map = Dict{Float64, Datum}()
    for (mu, arr) in pairs(unaveraged_data)
        total_sign = weightedmean([datum.total_sign for datum in arr])
        avg_energy = weightedmean([datum.avg_energy for datum in arr])
        n = weightedmean([datum.n for datum in arr])
        nup = weightedmean([datum.nup for datum in arr])
        ndn = weightedmean([datum.ndn for datum in arr])
        nupndn = weightedmean([datum.nupndn for datum in arr])
        z_spin = (nup - ndn) / 2
        moment = n - 2 * nupndn
        ferro_cf = weightedmean([datum.ferro_cf for datum in arr])
        antiferro_cf = weightedmean([datum.antiferro_cf for datum in arr])

        connected = [n^2, n * nupndn, nupndn^2, moment^2, z_spin^2]

        cfs = Array{Matrix{Measurement}, 1}(undef, 5)
        for i in 1:5
            width = 2*(size(arr[1].cfs[1], 1)-1)
            cfs[i] = Matrix{Measurement}(undef, width, width)
            for j in 1:size(arr[1].cfs[1], 1)
                for k in j:size(arr[1].cfs[1], 1)
                    cfs[i][k, j] = weightedmean([datum.cfs[i][k, j] for datum in arr]) - connected[i]
                    cfs[i][j, k] = cfs[i][k, j]
                end
            end
        end
        datum_map[mu] = Datum(total_sign, avg_energy, n, nup, ndn, nupndn, z_spin, moment,
                              ferro_cf, antiferro_cf, cfs)
    end
    return datum_map
end

function flesh_out_cfs!(averaged_data)
    for (mu, datum) in averaged_data
        n = Int(size(datum.cfs[1], 1)/2)
        fulln = size(datum.cfs[1], 1)
        for i in 1:5
            for x in 0:n
                for y in 0:n
                    datum.cfs[i][mod(-y, fulln) + 1, x + 1] =
                        datum.cfs[i][y + 1, x + 1]
                    datum.cfs[i][mod(-y, fulln) + 1, mod(-x, fulln) + 1] =
                        datum.cfs[i][y + 1, x + 1]
                    datum.cfs[i][y + 1, mod(-x, fulln) + 1] =
                        datum.cfs[i][y + 1, x + 1]
                end
            end
        end
    end
end

function calculate_SAFs(averaged_data)
    # (pi,pi) fourier-transformed spin correlation function
    SAFs = Dict{Float64, Measurement{Float64}}()
    qx = pi
    qy = pi
    for (mu, datum) in averaged_data
        SAFs[mu] = 0
        spin_cf = datum.cfs[5]
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

function run()
    for beta_dir in filter(isdir, readdir("10x10"; join = true))
        println(beta_dir)
        cd(beta_dir)
        unaveraged_data = Dict{Float64, Array{Datum, 1}}()
        beta = 0
        for filename in readdir()
            if filename[1:2] != "rz"
                continue
            end
            header = parse_header(filename)
            beta = header.beta
            #println(header)
            datum = parse_file(filename)
            get!(unaveraged_data, header.mu, Array{Datum, 1}(undef, 0))
            push!(unaveraged_data[header.mu], datum)
        end
        datum_map = average_data(unaveraged_data)
        flesh_out_cfs!(datum_map)
        safs = calculate_SAFs(datum_map)
        fname = "../../safs_beta_$beta.dat"
        file = open(fname, "w")
        for (mu, saf) in safs
            n = Measurements.value(datum_map[mu].n)
            n_uc = Measurements.uncertainty(datum_map[mu].n)
            saf = Measurements.value(saf)
            saf_uc = Measurements.uncertainty(saf)
            write(file, "$n $n_uc  $saf  $saf_uc\n")
        end
        close(file)
        cd("../../")
    end
end

run()
