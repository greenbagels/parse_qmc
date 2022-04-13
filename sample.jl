
struct Header
    n::Int
    l::Int
    u::Float64
    dt::Float64
    mu::Float64
    beta::Float64
end

# Size: 10 + 2 * (5 * n^2) = 10(n^2 + 1) * sizeof(float64), doubled to count
#       uncertainties.

# So, for example, when n=10, this is just 16160 bytes
# If we tag with an additional mu index, it's 16168 bytes
struct Sample{T}
    total_sign::Measurement{T}
    avg_energy::Measurement{T}
    n::Measurement{T}
    nup::Measurement{T}
    ndn::Measurement{T}
    nupndn::Measurement{T}
    z_spin::Measurement{T}
    moment::Measurement{T}
    ferro_cf::Measurement{T}
    antiferro_cf::Measurement{T}
    cfs::Array{Matrix{Measurement{T}}, 1}
    sfs::Array{Matrix{Measurement{T}}, 1}
end
