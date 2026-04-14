# AB stacked bilayer graphene BC and Chern number calculation

using Plots
using LinearAlgebra
using Measures # for nicer plot padding

# hopping amplitudes
t1 = 1; # intralayer hopping
t2 = 0.2; # interlayer hopping
Δ = 0.5; # sublattice potential difference

# lattice vectors
a1 = (sqrt(3)/2).*[sqrt(3),1]
a2 = (sqrt(3)/2).*[sqrt(3),-1]

# bloch hamiltonian
function get_h(kx::Real, ky::Real, t1::Real = t1, t2::Real = t2, a1::Vector{Float64} = a1, a2::Vector{Float64} = a2, Δ::Real = Δ)
    f1 = t1*(1 + exp(-im*(a1[1]*kx + a1[2]*ky)) + exp(-im*(a2[1]*kx + a2[2]*ky)))
    f1c = conj(f1)
    h = [Δ 0 f1c 0;
         0 Δ 0 f1;
         f1 0 -Δ t2;
         0 f1c t2 -Δ]
    return Hermitian(h)
end

# -------------------- plot the band structure
# path in BZ: K (hexagon corner) → Γ (origin) → M (middle of BZ edge) → K' (other corner)
Nkpath = 100 # number of points in each path section
kpath = [(0.0, (1-j)*4π/3/sqrt(3)) for j in range(0,1,Nkpath)]
append!(kpath, [(2π/3 *j, 0.0) for j in range(0,1,Nkpath)])
append!(kpath, [(2π/3, j* 2π/3/sqrt(3)) for j in range(0,1,Nkpath)])

# calculation along path
energies = Array{Float64}(undef, Nkpath*3, 4)
for (nk,(kx,ky)) in enumerate(kpath)
    eigenenergies = eigvals(get_h(kx,ky))
    for (j,ee) in enumerate(eigenenergies)
        energies[nk,j] = ee
    end
end
bandplot = plot(framestyle=:box, ylabel="Energy", xlabel="", title = "Band Structure; t = $t1, t' = $t2",
            xticks=([0,Nkpath,2*Nkpath,3*Nkpath],["K", "Γ", "M", "K'"]), legend=:topright);
for j in 1:4
    en = energies[:,j]
    plot!(bandplot, en, label="");
end
# --------------------

# -------------------- calculate and plot the Berry curvature
# derivatives of hamiltonian
function get_dxh(kx::Real, ky::Real, t1::Real = t1, t2::Real = t2, a1::Vector{Float64} = a1, a2::Vector{Float64} = a2)
    e1 = exp(-im * (a1[1]*kx + a1[2]*ky))
    e2 = exp(-im * (a2[1]*kx + a2[2]*ky))
    df1 = t1 * (-im*a1[1]*e1 - im*a2[1]*e2)
    df1c = conj(df1)
    dh = [0   0   df1c  0;
          0   0   0     df1;
          df1 0   0     0;
          0   df1c 0    0]
    return Hermitian(dh)
end
function get_dyh(kx::Real, ky::Real, t1::Real = t1, t2::Real = t2, a1::Vector{Float64} = a1, a2::Vector{Float64} = a2)
    e1 = exp(-im * (a1[1]*kx + a1[2]*ky))
    e2 = exp(-im * (a2[1]*kx + a2[2]*ky))
    df1 = t1 * (-im*a1[2]*e1 - im*a2[2]*e2)
    df1c = conj(df1)
    dh = [0   0   df1c  0;
          0   0   0     df1;
          df1 0   0     0;
          0   df1c 0    0]
    return Hermitian(dh)
end

# Berry curvature of band n at kx, ky
function get_F(n::Integer, kx::Real,ky::Real)
    eigvals, eigvecs = eigen(get_h(kx,ky))
    F = 0.0
    for m in 1:4
        if m == n
            continue
        end
        Fm = eigvecs[:,n]' * get_dxh(kx,ky) * eigvecs[:,m]
        Fm *= eigvecs[:,m]' * get_dyh(kx,ky) * eigvecs[:,n]
        Fm /= (eigvals[n] - eigvals[m])^2
        Fm = imag(Fm)
        Fm *= -2
        F += Fm
    end
    return F
end

# discretize the BZ
Nk = 500 # number of steps in kx and ky each
kx_range = range(-4π/3/sqrt(3),4π/3/sqrt(3),Nk)
ky_range = copy(kx_range)
dk = step(kx_range)

F_array = Array{Float64}(undef,Nk,Nk)
for (i,kx) in enumerate(kx_range)
    for (j,ky) in enumerate(ky_range)
        F_array[i,j] = get_F(1,kx,ky) # lowest band
    end
end

# calculate the Chern number
C = round(sum(F_array)*dk^2 / (2π); digits = 2)

bcplot = heatmap(kx_range,ky_range,F_array', title = "Ωᶻ → C = $C");
# ---------------------


# output everything in one plot
plt = plot(bandplot,bcplot,layout = (1,2), size = (1000,400), margins = 5mm, dpi=200);
display(plt)
savefig(plt, "bilayer_graphene_Delta$Δ.png")
