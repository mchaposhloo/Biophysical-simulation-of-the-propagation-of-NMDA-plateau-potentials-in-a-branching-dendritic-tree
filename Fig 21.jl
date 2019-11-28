
using DifferentialEquations
using Plots
using BlockArrays
using SparseArrays
using IterativeSolvers

GR = zeros(1,39)
amplitudes = zeros(3,39);

for counter in 1:39    
    r_m = 50 # kOhm cm^2 - specific membrane resistance (10000 Ωcm^2 Major et al, 2008) 
    c_m = 1 # μFarad/cm^2 - specific membrane capacitance (c) (0.8 Major et al, 2008)
    r_L = 0.1 # kOhm cm - longitudinal resistivity
    L = 0.04 #cm - length of each branch
    Δx = 0.001 # cm
    Δt = 0.01 #ms
    SOMA_SURFACE_AREA = 0.00002 # cm^2
    C_s = SOMA_SURFACE_AREA*c_m # μFarad - somatic membrane capacitance
    τ_m = r_m*c_m; # ms - membrane time constant

    d_0_apical_trunk = 0.0008 #cm - diameter at the begining of the apical trunk
    #d_avg_apical_trunk = 0.0006 #cm - average diameter across the apical trunk
    d_L_apical_trunk = 0.0006; #cm - diameter at the end of the apical trunk
    #defining the linear diameter function for the apical trunk
    d(x) = ((d_L_apical_trunk - d_0_apical_trunk)/L)*x + d_0_apical_trunk;

    d_0_branch_1 = 0.0001 #cm - diameter at the begining of the 1st branch
    d_L_branch_1 = 0.00005 #cm - diameter at the end of the 1st branch
    #the linear diameter function for the 1st branch
    d2(x) = ((d_L_branch_1 - d_0_branch_1)/L)*x + d_0_branch_1;

    d_0_branch_2 = 0.0006 #cm - diameter at the begining of the 2st branch
    d_L_branch_2 = 0.0004 #cm - diameter at the end of the 2st branch
    #the linear diameter function for the 2st branch
    d3(x) = ((d_L_branch_2 - d_0_branch_2)/L)*x + d_0_branch_2;

    d_0_branch_3_4 = 0.0001 #cm - diameter at the begining of the 2st branch
    d_L_branch_3_4 = 0.00005 #cm - diameter at the end of the 2st branch
    #the linear diameter function for the 2st branch
    d4(x) = ((d_L_branch_3_4 - d_0_branch_3_4)/L)*x + d_0_branch_3_4;
    d_L = 0.00001 + (counter-1)*0.00001
    d_0 = 0.0004 - (counter-1)*0.00001
    d5(x) = ((d_L - d_0)/L)*x + d_0;
    GR[counter] = (0.0004^(3/2) + 0.00005^(3/2))/((d_0)^(3/2))

    shared_node = Int(round(L/Δx))+1
    terminal_node_1 = shared_node + shared_node - 1
    shared_node_2 = terminal_node_1 + shared_node - 1
    terminal_node_3 = shared_node_2 + shared_node - 1
    terminal_node_4 = terminal_node_3 + shared_node - 1

    A = BlockArray(zeros(terminal_node_4, terminal_node_4), [1, shared_node - 2, 1, terminal_node_1 - shared_node, shared_node_2 - terminal_node_1 - 1, 1, terminal_node_3 - shared_node_2, terminal_node_4 - terminal_node_3], [1, shared_node - 2, 1, terminal_node_1 - shared_node, shared_node_2 - terminal_node_1 - 1, 1, terminal_node_3 - shared_node_2, terminal_node_4 - terminal_node_3]);

    block_parent = zeros(shared_node-2,shared_node-2);
    α = Δt/(4*r_L*c_m*(Δx^2));
    d_square_half_minus(i,d_func) = ((d_func(i*Δx))^2 + (d_func((i-1)*Δx))^2)/2
    d_square_half_plus(i,d_func) = ((d_func(i*Δx))^2 + (d_func((i+1)*Δx))^2)/2
    a(i,d_func) = -α*d_square_half_minus(i,d_func)/d_func(i*Δx)
    b(i,d_func) = 1 + (Δt/τ_m) + (α*(d_square_half_minus(i,d_func) + d_square_half_plus(i,d_func))/d_func(i*Δx))
    c(i,d_func) = -α*d_square_half_plus(i,d_func)/d_func(i*Δx);

    for i in 1:size(block_parent)[1] #excluding the soma and the boundary for now
        for j in 1:size(block_parent)[2]
            if (abs(j-i) <= 1) && (i > j)
                block_parent[i,j] = a(i,d)
            elseif (abs(j-i) <= 1) && (i < j)
                block_parent[i,j] = c(i,d)
            elseif i == j
                block_parent[i,j] = b(i,d)
            end
        end
    end

    setblock!(A, block_parent, 2, 2);

    block_child_1 = zeros(terminal_node_1 - shared_node,terminal_node_1 - shared_node);

    for i in 1:size(block_child_1)[1] 
        for j in 1:size(block_child_1)[2]
            if (i == size(block_child_1)[1]) && (j == (size(block_child_1)[2]-1))
                block_child_1[i,j] = a(i,d2) + c(i,d2)
            elseif (i == size(block_child_1)[1]) && (j == size(block_child_1)[2])
                block_child_1[i,j] = b(i,d2)
            elseif (abs(j-i) <= 1) && (i > j)
                block_child_1[i,j] = a(i,d2)
            elseif (abs(j-i) <= 1) && (i < j)
                block_child_1[i,j] = c(i,d2)
            elseif i == j
                block_child_1[i,j] = b(i,d2)
            end
        end
    end

    setblock!(A, block_child_1, 4, 4);

    block_child_2 = zeros(shared_node_2 - terminal_node_1 - 1,shared_node_2 - terminal_node_1 - 1);

    for i in 1:size(block_child_2)[1] #excluding the soma and the boundary for now
        for j in 1:size(block_child_2)[2]
            if (abs(j-i) <= 1) && (i > j)
                block_child_2[i,j] = a(i,d3)
            elseif (abs(j-i) <= 1) && (i < j)
                block_child_2[i,j] = c(i,d3)
            elseif i == j
                block_child_2[i,j] = b(i,d3)
            end
        end
    end

    setblock!(A, block_child_2, 5, 5);

    block_child_3 = zeros(terminal_node_3 - shared_node_2, terminal_node_3 - shared_node_2)

    for i in 1:size(block_child_3)[1] 
        for j in 1:size(block_child_3)[2]
            if (i == size(block_child_3)[1]) && (j == (size(block_child_3)[2]-1))
                block_child_3[i,j] = a(i,d4) + c(i,d4)
            elseif (i == size(block_child_3)[1]) && (j == size(block_child_3)[2])
                block_child_3[i,j] = b(i,d4)
            elseif (abs(j-i) <= 1) && (i > j)
                block_child_3[i,j] = a(i,d4)
            elseif (abs(j-i) <= 1) && (i < j)
                block_child_3[i,j] = c(i,d4)
            elseif i == j
                block_child_3[i,j] = b(i,d4)
            end
        end
    end

    setblock!(A, block_child_3, 7, 7);

    block_child_4 = zeros(terminal_node_4 - terminal_node_3,terminal_node_4 - terminal_node_3)

    for i in 1:size(block_child_4)[1] #excluding the soma and the boundary for now
        for j in 1:size(block_child_4)[2]
            if (i == size(block_child_4)[1]) && (j == (size(block_child_4)[2]-1))
                block_child_4[i,j] = a(i,d5) + c(i,d5)
            elseif (i == size(block_child_2)[1]) && (j == size(block_child_2)[2])
                block_child_4[i,j] = b(i,d5)
            elseif (abs(j-i) <= 1) && (i > j)
                block_child_4[i,j] = a(i,d5)
            elseif (abs(j-i) <= 1) && (i < j)
                block_child_4[i,j] = c(i,d5)
            elseif i == j
                block_child_4[i,j] = b(i,d5)
            end
        end
    end

    setblock!(A, block_child_4, 8, 8);


    A = Array(A);

    β(i,d_func) = 2*r_L*c_m*d_func(i*Δx)*Δx
    λ(i,d_func) = (-C_s/Δt - π*(d_func(i*Δx))^2/(4*r_L)*β(i,d_func)/(Δt*d_square_half_minus(i,d_func)))^-1

    A[1,1] = λ(0,d)*(-C_s/Δt - SOMA_SURFACE_AREA/r_m
        + π*(d(0))^2/(4*r_L)*(-β(0,d)/(Δt*d_square_half_minus(0,d))
            - β(0,d)/(τ_m*d_square_half_minus(0,d)) 
            - 1/(2*Δx)*(1 + (d_square_half_plus(0,d) / d_square_half_minus(0,d)))))

    A[1,2] = λ(0,d)*π*(d(0))^2/(4*r_L)/(2*Δx)*(1 + (d_square_half_plus(0,d) / d_square_half_minus(0,d)))

    A[2,1] = a(1,d)

    γ = (-(d(shared_node*Δx))^2*β(shared_node,d)/Δt/d_square_half_plus(shared_node,d)
        - (d2(0))^2*β(0,d2)/Δt/d_square_half_minus(0,d2) 
        - (d3(0))^2*β(0,d3)/Δt/d_square_half_minus(0,d3))^-1

    A[shared_node,shared_node] = γ*(-(d(shared_node*Δx))^2*((β(shared_node,d)/(Δt*d_square_half_plus(shared_node,d)) + β(shared_node,d)/(τ_m*d_square_half_plus(shared_node,d)) + 1/(2*Δx)*(1 + (d_square_half_minus(shared_node,d) / d_square_half_plus(shared_node,d)))))) + γ*(d2(0))^2*(-(β(0,d2)/(Δt*d_square_half_minus(0,d2))) - β(0,d2)/(τ_m*d_square_half_minus(0,d2)) - 1/(2*Δx)*(1 + (d_square_half_plus(0,d2) / d_square_half_minus(0,d2)))) + γ*(d3(0))^2*(-(β(0,d3)/(Δt*d_square_half_minus(0,d3))) - β(0,d3)/(τ_m*d_square_half_minus(0,d3)) - 1/(2*Δx)*(1 + (d_square_half_plus(0,d3) / d_square_half_minus(0,d3))))

    A[shared_node,shared_node-1] = γ*(d(shared_node*Δx))^2/(2*Δx)*(d_square_half_minus(shared_node,d)/d_square_half_plus(shared_node,d)+1)

    A[shared_node-1,shared_node] = c(shared_node-1,d)

    A[shared_node,shared_node+1] = γ*(d2(0))^2/(2*Δx)*(d_square_half_plus(0,d2)/d_square_half_minus(0,d2) + 1)

    A[shared_node+1,shared_node] = a(1,d2)

    A[shared_node,terminal_node_1+1] = γ*(d3(0))^2/(2*Δx)*(d_square_half_plus(0,d3)/d_square_half_minus(0,d3) + 1)

    A[terminal_node_1+1,shared_node] = a(1,d3);


    γ_2 = (-(d3((shared_node_2-terminal_node_1)*Δx))^2*β((shared_node_2-terminal_node_1),d3)/Δt/d_square_half_plus((shared_node_2-terminal_node_1),d3)
        - (d4(0))^2*β(0,d4)/Δt/d_square_half_minus(0,d4) 
        - (d5(0))^2*β(0,d5)/Δt/d_square_half_minus(0,d5))^-1

    A[shared_node_2,shared_node_2] = γ_2*(-(d3((shared_node_2-terminal_node_1)*Δx))^2*((β((shared_node_2-terminal_node_1),d3)/(Δt*d_square_half_plus((shared_node_2-terminal_node_1),d3)) + β((shared_node_2-terminal_node_1),d3)/(τ_m*d_square_half_plus((shared_node_2-terminal_node_1),d3)) + 1/(2*Δx)*(1 + (d_square_half_minus((shared_node_2-terminal_node_1),d3) / d_square_half_plus((shared_node_2-terminal_node_1),d3)))))) + γ_2*(d4(0))^2*(-(β(0,d4)/(Δt*d_square_half_minus(0,d4))) - β(0,d4)/(τ_m*d_square_half_minus(0,d4)) - 1/(2*Δx)*(1 + (d_square_half_plus(0,d4) / d_square_half_minus(0,d4)))) + γ_2*(d5(0))^2*(-(β(0,d5)/(Δt*d_square_half_minus(0,d5))) - β(0,d5)/(τ_m*d_square_half_minus(0,d5)) - 1/(2*Δx)*(1 + (d_square_half_plus(0,d5) / d_square_half_minus(0,d5))))

    A[shared_node_2,shared_node_2-1] = γ_2*(d3((shared_node_2-terminal_node_1)*Δx))^2/(2*Δx)*(d_square_half_minus((shared_node_2-terminal_node_1),d3)/d_square_half_plus((shared_node_2-terminal_node_1),d3)+1)

    A[shared_node_2-1,shared_node_2] = c(shared_node_2-terminal_node_1-1,d3)

    A[shared_node_2,shared_node_2+1] = γ_2*(d4(0))^2/(2*Δx)*(d_square_half_plus(0,d4)/d_square_half_minus(0,d4) + 1)

    A[shared_node_2+1,shared_node_2] =a(1,d4)

    A[shared_node_2,terminal_node_3+1] = γ_2*(d5(0))^2/(2*Δx)*(d_square_half_plus(0,d5)/d_square_half_minus(0,d5) + 1)

    A[terminal_node_3+1,shared_node_2] = a(1,d5);

    A_sparse = sparse(A);

    function solve_linear(A,U,clamp_pos_index, vol_clamp_val, onset, dur)
        for t in 1:(size(U,2)-1)
            U[:,t+1] .= gmres(A,U[:,t],tol = 1e-5) #A\U[:,t]
            if (t>=(round(onset/Δt)+1)) && (t<=(round((onset+dur)/Δt)+1))
                U[clamp_pos_index,t+1] = vol_clamp_val
            end
        end
        return U
    end
    M_t = 15000;
    V = zeros(terminal_node_4,M_t+1);

    V_clamp_4th_cild = solve_linear(A_sparse,V,162,40,0.01,100);
    amplitudes[1,counter] = maximum(V_clamp_4th_cild[162,:])
    amplitudes[2,counter] = maximum(V_clamp_4th_cild[122,:])
    amplitudes[3,counter] = maximum(V_clamp_4th_cild[120,:])
end

plot(log.(GR)[:],amplitudes[1,:],seriestype=:scatter,markershape = :hexagon,markersize = 2.5,label = "3rd child")
plot!(log.(GR)[:],amplitudes[2,:],seriestype=:scatter,markershape = :square,markersize = 2.5,label = "4th child")
plot!(log.(GR)[:],amplitudes[3,:],seriestype=:scatter,markershape = :circle,markersize = 2.5,label = "2nd child",xaxis ="log(GR)",yaxis ="voltage (mv)",legend = :bottomleft,ylims = (0,41))#ylims = (15.5,35.5),, legendfontsize =8

savefig("GR voltage drop 1.svg")

plot(log.(GR)[:],((amplitudes[1,:]-amplitudes[2,:])./amplitudes[1,:])*100,seriestype=:scatter,markershape = :hexagon,label = "",xaxis ="log(GR)",yaxis ="% voltage drop",ylims = (0,86),markersize = 2.5)

savefig("GR voltage drop 2.svg")
