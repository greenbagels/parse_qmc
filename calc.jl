using BenchmarkTools
using Plots
using Plots.PlotMeasures
using FiniteDifferences

include("parse.jl")

function parse_beta_directory()
    unaveraged_data = Dict{Float64, Array{Sample{Float64}, 1}}()
    beta = 0
    for filename in readdir()
        if filename[1:2] != "rz"
            continue
        end
        header = parse_bandmott_header(filename)
        beta = header.beta
        sample = parse_bandmott_file(filename)
        get!(unaveraged_data, header.mu, Array{Sample{Float64}, 1}(undef, 0))
        push!(unaveraged_data[header.mu], sample)
    end
    return (beta, unaveraged_data)
end

function print_faf_sf_dat(beta_str, samples)
    println("Printing SF (Antiferro) to file")

    fname = "antiferro_beta_$beta_str.dat"
    open(fname, "w") do file
        for (mu, sample) in sort(collect(samples), by=x->x[1])
            println("starting loop for mu = $mu")
            n = Measurements.value(sample.n)
            n_uc = Measurements.uncertainty(sample.n)
            temp = calculate_single_sf(sample.cfs[5], pi, pi)
            # saf = Measurements.value(sample.sfs[5][5, 5])
            saf = Measurements.value(temp)
            # saf_uc = Measurements.uncertainty(sample.sfs[5][5, 5])
            saf_uc = Measurements.uncertainty(temp)
            write(file, "$n $n_uc  $saf  $saf_uc\n")
        end
    end

    println("Printing SF (Ferro) to file")

    fname = "ferro_beta_$beta_str.dat"
    open(fname, "w") do file
        for (mu, sample) in sort(collect(samples), by=x->x[1])
            println("starting loop for mu = $mu")
            n = Measurements.value(sample.n)
            n_uc = Measurements.uncertainty(sample.n)
            temp = calculate_single_sf(sample.cfs[5], 0, 0)
            # saf = Measurements.value(sample.sfs[5][5, 5])
            saf = Measurements.value(temp)
            # saf_uc = Measurements.uncertainty(sample.sfs[5][5, 5])
            saf_uc = Measurements.uncertainty(temp)
            write(file, "$n $n_uc  $saf  $saf_uc\n")
        end
    end
end

function print_sf_heatmap(beta_str, samples)
    mkpath("sf_heatmap_data/$beta_str/img/spin")
    mkpath("sf_heatmap_data/$beta_str/img/moment")
    mkpath("sf_heatmap_data/$beta_str/dat/spin")
    mkpath("sf_heatmap_data/$beta_str/dat/moment")
    cd("sf_heatmap_data/$beta_str")
    i = 1
    for (mu, sample) in sort(collect(samples), by=x->x[1])
        index = lpad(i, 4, "0")
        i += 1

        open("dat/moment/sf_$index.dat", "w") do file
            write(file, "# u = 8, l = 10, n = $(Measurements.value(sample.n))\n")
            for iy in 1:10
                for ix in 1:10
                    z = sample.sfs[4][ix, iy]
                    write(file, "$ix $iy $(Measurements.value(z))\n")
                end
                write(file, "\n");
            end
        end

        open("dat/spin/sf_$index.dat", "w") do file
            write(file, "# u = 8, l = 10, n = $(Measurements.value(sample.n))\n")
            for iy in 1:10
                for ix in 1:10
                    z = sample.sfs[5][ix, iy]
                    write(file, "$ix $iy $(Measurements.value(z))\n")
                end
                write(file, "\n");
            end
        end

        qs = [x * 2. * pi / 10. for x in 1:10]
        n = round(Measurements.value(sample.n), digits=4)
        heatmap(qs, qs, map(Measurements.value, sample.sfs[5]),
                xticks = 0:2*pi:pi/4, yticks = 0:2*pi:pi/4, #clims=(0, 2),
                # title = "\$S_z Structure Factor, l=10, U=8, \\beta=$beta_str, n=$n\$",
                xlabel = "q_x", ylabel = "q_y",
                left_margin = 15px, right_margin = 15px, bottom_margin = 15px)
        savefig("img/spin/sf_$index.png")
        heatmap(qs, qs, map(Measurements.value, sample.sfs[4]), #clims=(0, 2),
                xticks = 0:2*pi:pi/4, yticks = 0:2*pi:pi/4,
                # title = "\$m_z^2 Structure Factor, l=10, U=8, \\beta=$beta_str, n=$n\$",
                xlabel = "q_x", ylabel = "q_y",
                left_margin = 15px, right_margin = 15px, bottom_margin = 15px)
        savefig("img/moment/sf_$index.png")
    end
    cd("../../")
end

function print_sf_dat(beta_str, samples)
    # unavged map: mu -> (r -> list of sf(r))
    spin_sf_map = Dict{Float64, Dict{Float64, Measurement{Float64}}}()
    moment_sf_map = Dict{Float64, Dict{Float64, Measurement{Float64}}}()
    for (mu, sample) in samples
        temp_spin_sf = Dict{Float64, Vector{Measurement{Float64}}}()
        temp_moment_sf = Dict{Float64, Vector{Measurement{Float64}}}()
        for i in 1:10
            for j in 1:10
                r = sqrt((i-1)^2 + (j-1)^2)
                get!(temp_spin_sf, r, Vector{Measurement{Float64}}(undef, 0))
                get!(temp_moment_sf, r, Vector{Measurement{Float64}}(undef, 0))
                push!(temp_spin_sf[r], sample.sfs[5][i, j])
                push!(temp_moment_sf[r], sample.sfs[4][i, j])
            end
        end
        spin_sf_map[mu] = Dict([(r, weightedmean(vec)) for (r, vec) in temp_spin_sf])
        moment_sf_map[mu] = Dict([(r, weightedmean(vec)) for (r, vec) in temp_moment_sf])
    end

    fname = "spin_sf_beta_$beta_str.dat"
    open(fname, "w") do file
        for (mu, rsf) in sort(collect(spin_sf_map), by=x->x[1])
            n = Measurements.value(samples[mu].n)
            n_uc = Measurements.uncertainty(samples[mu].n)
            write(file, "$mu $n $n_uc ")
            for (r, sf) in sort(collect(rsf), by=x->x[1])
                val = Measurements.value(sf)
                uc = Measurements.uncertainty(sf)
                write(file, "$val $uc ")
            end
            write(file, "\n")
        end
    end
    fname = "moment_sf_beta_$beta_str.dat"
    open(fname, "w") do file
        for (mu, rsf) in sort(collect(moment_sf_map), by=x->x[1])
            n = Measurements.value(samples[mu].n)
            n_uc = Measurements.uncertainty(samples[mu].n)
            write(file, "$mu $n $n_uc ")
            for (r, sf) in sort(collect(rsf), by=x->x[1])
                val = Measurements.value(sf)
                uc = Measurements.uncertainty(sf)
                write(file, "$val $uc ")
            end
            write(file, "\n")
        end
    end
end

function run()
    ENV["GKSwstype"] = "nul"
    gr()
    benchmark = false
    reuse_old = true
    samples = Dict{Float64, Dict{Float64, Sample{Float64}}}()
    for beta_dir in filter(isdir, readdir("10x10"; join = true))

        beta_str = match(r"beta([0-9]+\.[0-9]+)", beta_dir).captures[1]
        beta = parse(Float64, beta_str)
        # println(beta_dir)
        cd(beta_dir)
        cd("../../")

        old_fname = "averaged_$beta_str.old"

        if isfile(old_fname) && reuse_old
            println("Found old averaging file: $old_fname ...")
            samples[beta] = read_avg_data_from_file(old_fname)
        else
            println("No old averaging file found! Parsing data files ourselves...")
            cd(beta_dir)

            if benchmark
            #    println("Benchmarking directory parsing...")
            #    @btime parse_beta_directory()
            end

            (beta, unaveraged_data) = parse_beta_directory()
            cd("../../")

            if benchmark
                println("Benchmarking data averaging...")
                @btime average_data($unaveraged_data)
            end

            samples[beta] = average_data(unaveraged_data)
            write_avg_data_to_file(old_fname, samples[beta])
        end

#       if benchmark
#           println("Benchmarking SAF calculation...")
#           @btime calculate_SAFs($sample_map)
#       end

        #print_faf_sf_dat(beta_str, samples[beta])
        #print_sf_heatmap(beta_str, samples[beta])
        #print_sf_dat(beta_str, samples[beta])
    end

    for beta in [1.0 2.0 3.0]
        step = 0.05
        open("beta_$(beta)_extras.dat", "w") do file
            write(file, "n    dE/dT    sign\n")
            barr = [beta + step*x for x in -2:2]
            for mu in keys(samples[beta])
                dEdT = -beta^2 * central_fdm(5, 1)((x) -> samples[x][mu].avg_energy, beta, step)
                write(file, "$(samples[beta][mu].n)     $dEdT    $(samples[beta][mu].total_sign)\n")
            end
        end

        step = 0.1
        open("beta_$(beta)_kappa.dat", "w") do file
            write(file, "n    dndmu    kappa    mu\n")
            for mu in -14.5:0.5:-10.5
                dndmu = central_fdm(3, 1)((x) -> samples[beta][round(x; digits=1)].n, mu, 5*step)
                n = samples[beta][mu].n
                kappa = dndmu / n^2
                write(file, "$n    $dndmu    $kappa    $mu\n")
            end
            dndmu = FiniteDifferenceMethod([-5, 0, 1], 1)((x) -> samples[beta][round(x; digits=1)].n, -10, step)
            n = samples[beta][-10].n
            kappa = dndmu / n^2
            write(file, "$n    $dndmu    $kappa    -10\n")
            for mu in -9.9:0.1:9.9
                dndmu = central_fdm(3, 1)((x) -> samples[beta][round(x; digits=1)].n, mu, step)
                n = samples[beta][mu].n
                kappa = dndmu / n^2
                write(file, "$n    $dndmu    $kappa    $mu\n")
            end
            dndmu = FiniteDifferenceMethod([-1, 0, 5], 1)((x) -> samples[beta][round(x; digits=1)].n, 10, step)
            n = samples[beta][10].n
            kappa = dndmu / n^2
            write(file, "$n    $dndmu    $kappa    10\n")
            for mu in 10.5:0.5:14.5
                dndmu = central_fdm(3, 1)((x) -> samples[beta][round(x; digits=1)].n, mu, 5*step)
                n = samples[beta][mu].n
                kappa = dndmu / n^2
                write(file, "$n    $dndmu    $kappa    $mu\n")
            end
        end
    end
end

run()
