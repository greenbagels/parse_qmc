
include("parse.jl")

function run()
    for beta_dir in filter(isdir, readdir("10x10"; join = true))
        println(beta_dir)
        cd(beta_dir)
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
        sample_map = average_data(unaveraged_data)
        flesh_out_cfs!(sample_map)
        safs = calculate_SAFs(sample_map)
        fname = "../../safs_beta_$beta.dat"
        file = open(fname, "w")
        for (mu, saf) in safs
            n = Measurements.value(sample_map[mu].n)
            n_uc = Measurements.uncertainty(sample_map[mu].n)
            saf = Measurements.value(saf)
            saf_uc = Measurements.uncertainty(saf)
            write(file, "$n $n_uc  $saf  $saf_uc\n")
        end
        close(file)
        cd("../../")
    end
end

run()
