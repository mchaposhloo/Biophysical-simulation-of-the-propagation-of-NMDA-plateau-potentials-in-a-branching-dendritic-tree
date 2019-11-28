
using DifferentialEquations
using Plots
using BlockArrays
using SparseArrays
using IterativeSolvers

r_m = 50 # kOhm cm^2 - specific membrane resistance (10000 Ωcm^2 Major et al, 2008) 
c_m = 1 # μFarad/cm^2 - specific membrane capacitance (c) (0.8 Major et al, 2008)
r_L = 0.1 # kOhm cm - longitudinal resistivity
L = 0.04 #cm - length of each branch
Δx = 0.001 # cm
Δt = 0.01 #ms
#SOMA_SURFACE_AREA = 0.00002 # cm^2
#C_s = SOMA_SURFACE_AREA*c_m # μFarad - somatic membrane capacitance
τ_m = r_m*c_m; # ms - membrane time constant



d_0_branch_3_4 = 0.0001 #cm - diameter at the begining of the 2st branch
d_L_branch_3_4 = 0.00005 #cm - diameter at the end of the 2st branch
#the linear diameter function for the 2st branch
#d4(x) = ((d_L_branch_3_4 - d_0_branch_3_4)/L)*x + d_0_branch_3_4;
d5(x) = ((d_L_branch_3_4 - d_0_branch_3_4)/L)*x + d_0_branch_3_4;

terminal_node = Int(round(L/Δx))+1



block = zeros(terminal_node,terminal_node);
α = Δt/(4*r_L*c_m*(Δx^2));
d_square_half_minus(i,d_func) = ((d_func(i*Δx))^2 + (d_func((i-1)*Δx))^2)/2
d_square_half_plus(i,d_func) = ((d_func(i*Δx))^2 + (d_func((i+1)*Δx))^2)/2
a(i,d_func) = -α*d_square_half_minus(i,d_func)/d_func(i*Δx)
b(i,d_func) = 1 + (Δt/τ_m) + (α*(d_square_half_minus(i,d_func) + d_square_half_plus(i,d_func))/d_func(i*Δx))
c(i,d_func) = -α*d_square_half_plus(i,d_func)/d_func(i*Δx);



for i in 1:size(block)[1] 
    for j in 1:size(block)[2]
        if i == 1
            block[i,i] = b(i,d5)
            block[i,i+1] = a(i,d5) + c(i,d5)
        elseif (i == size(block)[1]) && (j == (size(block)[2]-1))
            block[i,j] = a(i,d5) + c(i,d5)
        elseif (i == size(block)[1]) && (j == size(block)[2])
            block[i,j] = b(i,d5)
        elseif (abs(j-i) <= 1) && (i > j)
            block[i,j] = a(i,d5)
        elseif (abs(j-i) <= 1) && (i < j)
            block[i,j] = c(i,d5)
        elseif i == j
            block[i,j] = b(i,d5)
        end
    end
end


A_sparse = sparse(block);

function solve_linear(A,U,clamp_pos_index, vol_clamp_val, onset, dur)
    for t in 1:(size(U,2)-1)
        U[:,t+1] .= gmres(A,U[:,t],tol = 1e-5) #A\U[:,t]
        if (t>=(round(onset/Δt)+1)) && (t<=(round((onset+dur)/Δt)+1))
            U[clamp_pos_index,t+1] = vol_clamp_val
        end
    end
    return U
end

M_t = 25000;
V = zeros(terminal_node,M_t+1);

V_clamp_4th_cild = solve_linear(A_sparse,V,20,40,0.01,100);

t = 0:Δt:(M_t*Δt);

p1 = plot(t,V_clamp_4th_cild[1,:],
    xaxis ="time (ms)",yaxis ="voltage (mv)", lw = 2.5,fmt = :svg,label = "beginning", legendfontsize =7)
p2 = plot(t,V_clamp_4th_cild[terminal_node,:],
    xaxis ="time (ms)", lw = 2.5,fmt = :svg,label = "end",linecolor = :red, legendfontsize =8)
plot(p1,p2,layout=(1,2))

savefig("sealed-end.svg")
