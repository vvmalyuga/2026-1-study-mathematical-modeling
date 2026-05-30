using DrWatson
@quickactivate "firm_competition"
using DifferentialEquations, Plots, DataFrames, JLD2

script_name = splitext(basename(@__FILE__))[1]
plots_dir = plotsdir(script_name)
data_dir = datadir(script_name)
mkpath(plots_dir)
mkpath(data_dir)

p_cr = 15.0     # критическая цена
N = 17.0        # число потребителей (тыс.)
q = 1.0
τ1 = 11.0
τ2 = 14.0
p̃1 = 8.0
p̃2 = 6.0
M1_0 = 2.5      # млн
M2_0 = 1.5

a1 = p_cr / (τ1^2 * p̃1^2 * N * q)
a2 = p_cr / (τ2^2 * p̃2^2 * N * q)
b = p_cr / (τ1^2 * p̃1^2 * τ2^2 * p̃2^2 * N * q)
c1 = (p_cr - p̃1) / (τ1 * p̃1)
c2 = (p_cr - p̃2) / (τ2 * p̃2)

println("Коэффициенты:")
println("a1 = $a1, a2 = $a2, b = $b, c1 = $c1, c2 = $c2")

function competition_case1!(du, u, p, t)
    M1, M2 = u
    du[1] = M1 - (b/c1)*M1*M2 - (a1/c1)*M1^2
    du[2] = (c2/c1)*M2 - (b/c1)*M1*M2 - (a2/c1)*M2^2
end

function competition_case2!(du, u, p, t)
    M1, M2 = u
    du[1] = M1 - (b/c1 + 0.001)*M1*M2 - (a1/c1)*M1^2
    du[2] = (c2/c1)*M2 - (b/c1)*M1*M2 - (a2/c1)*M2^2
end

θ_span = (0.0, 30.0)      # безразмерное время θ = t * c1
u0 = [M1_0, M2_0]

prob1 = ODEProblem(competition_case1!, u0, θ_span)
sol1 = solve(prob1, Tsit5(); saveat=0.1)

prob2 = ODEProblem(competition_case2!, u0, θ_span)
sol2 = solve(prob2, Tsit5(); saveat=0.1)

p1 = plot(sol1, idxs=[1,2],
          label=["Фирма 1" "Фирма 2"],
          xlabel="θ (безразмерное время)", ylabel="Оборотные средства M (млн)",
          title="Случай 1: только рыночная конкуренция",
          lw=2, legend=:bottomright)
savefig(p1, joinpath(plots_dir, "case1_curves.png"))

p2 = plot(sol2, idxs=[1,2],
          label=["Фирма 1" "Фирма 2"],
          xlabel="θ", ylabel="M (млн)",
          title="Случай 2: с социально-психологическим фактором (0.001)",
          lw=2, legend=:bottomright)
savefig(p2, joinpath(plots_dir, "case2_curves.png"))

A = [a1/c1   b/c1;
     b/c1    a2/c1]
b_vec = [1.0, c2/c1]

M_stat = A \ b_vec   # решение линейной системы
M1_stat, M2_stat = M_stat[1], M_stat[2]

println("\nСтационарное состояние (случай 1):")
println("M1* = ", M1_stat)
println("M2* = ", M2_stat)

@save datadir(script_name, "solution_case1.jld2") sol1
@save datadir(script_name, "solution_case2.jld2") sol2
@save datadir(script_name, "stationary_case1.jld2") M1_stat M2_stat

println("\n✅ Работа завершена. Графики сохранены в $plots_dir")
