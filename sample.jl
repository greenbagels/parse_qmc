
struct Header
    n::Int
    l::Int
    u::Float64
    dt::Float64
    mu::Float64
    beta::Float64
end

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
end
