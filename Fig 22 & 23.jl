
using DifferentialEquations
using Plots

function new_fitzhugh_nagumo_40!(du,u,p,t)
    a,b,c = p
    du[1] = u[1]*(a-u[1]/40)*(u[1]/40-1) - u[2]
    du[2] = b*u[1] - c*u[2]
end

u0 = [9,0.0]
p = (0.2, 0.00204,0.01)
tspan = (0.0, 400.0)
prob = ODEProblem(new_fitzhugh_nagumo_40!, u0, tspan, p)
sol1 = solve(prob)

u0 = [9,0.0]
p = (0.2, 0.00066,0.0039)
tspan = (0.0, 400.0)
prob = ODEProblem(new_fitzhugh_nagumo_40!, u0, tspan, p)
sol2 = solve(prob)

u0 = [9,0.0]
p = (0.2, 0.0016001,0.01)
tspan = (0.0, 400.0)
prob = ODEProblem(new_fitzhugh_nagumo_40!, u0, tspan, p)
sol3 = solve(prob)

plot(sol1,vars=(1), lw = 2,fmt = :svg,label = "b = 0.00204")
plot!(sol2,vars=(1),xaxis ="time (ms)",yaxis ="voltage (mv)", lw = 2,fmt = :svg,label = "b = 0.00175")
plot!(sol3,vars=(1), lw = 2,fmt = :svg,label = "b = 0.00160",xaxis ="time (ms)",yaxis ="voltage (mv)")

savefig("plateau FHN durations.svg")

u0 = [8.55,0.0]
p = (0.2, 0.00204,0.01)
tspan = (0.0, 200.0)
prob = ODEProblem(new_fitzhugh_nagumo_40!, u0, tspan, p)
sol2 = solve(prob)

u0 = [8,0.0]
p = (0.2, 0.00204,0.01)
tspan = (0.0, 200.0)
prob = ODEProblem(new_fitzhugh_nagumo_40!, u0, tspan, p)
sol3 = solve(prob)

plot(sol1,vars=(1), lw = 2,fmt = :svg,label = "plateau")
plot!(sol2,vars=(1),xaxis ="time (ms)",yaxis ="voltage (mv)", lw = 2,fmt = :svg,label = "subthreshold response")
plot!(sol3,vars=(1), lw = 2,fmt = :svg,label = "no response",xaxis ="time (ms)",yaxis ="voltage (mv)")

savefig("plateau FHN.svg")
