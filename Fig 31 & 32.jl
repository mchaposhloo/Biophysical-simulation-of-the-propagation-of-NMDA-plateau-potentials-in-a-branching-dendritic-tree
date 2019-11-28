
using DifferentialEquations
using Plots
using BlockArrays
using SparseArrays
using IterativeSolvers

r_m = 50 # kOhm cm^2 - specific membrane resistance (10000 Ωcm^2 Major et al, 2008) 
c_m = 1 # μFarad/cm^2 - specific membrane capacitance (c) (0.8 Major et al, 2008)
r_L = 0.1 # kOhm cm - longitudinal resistivity
L = 0.04 # cm - length of the cable
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

d_0_branch_1 = 0.0004 #cm - diameter at the begining of the 1st branch
d_L_branch_1 = 0.0002 #cm - diameter at the end of the 1st branch
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
d5(x) = ((d_L_branch_3_4 - d_0_branch_3_4)/L)*x + d_0_branch_3_4;

shared_node = Int(round(L/Δx))+1
terminal_node_1 = shared_node + shared_node - 1
shared_node_2 = terminal_node_1 + shared_node - 1
terminal_node_3 = shared_node_2 + shared_node - 1
terminal_node_4 = terminal_node_3 + shared_node - 1

A = BlockArray(zeros(terminal_node_4+2, terminal_node_4+2), [1, shared_node - 2, 1, terminal_node_1 - shared_node, shared_node_2 - terminal_node_1 - 1, 1, terminal_node_3 - shared_node_2, terminal_node_4 - terminal_node_3,2], [1, shared_node - 2, 1, terminal_node_1 - shared_node, shared_node_2 - terminal_node_1 - 1, 1, terminal_node_3 - shared_node_2, terminal_node_4 - terminal_node_3,2]);

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

A[terminal_node_4+1,terminal_node_4+1] = 0.2*Δt + 1
A[terminal_node_4+1,terminal_node_4+2] = Δt
A[terminal_node_4+2,terminal_node_4+1] = -0.00204*Δt
A[terminal_node_4+2,terminal_node_4+2] = 0.01*Δt + 1;

R_spine_mem = 200;

Spine_node = terminal_node_4-1;
Original_value = A[Spine_node,Spine_node]
A[Spine_node,Spine_node] = A[Spine_node,Spine_node] + Δt/(pi*c_m*d5(L)*R_spine_mem);
A[Spine_node,terminal_node_4+1] = -Δt/(pi*c_m*d5(L)*R_spine_mem);

A_sparse = sparse(A);

f(V_Spine_previous) = -Δt/(40^2)*V_Spine_previous^3 + Δt*(0.2+1)/40*V_Spine_previous^2;



function solve_linear(A,U,I)
    for t in 1:(size(U,2)-1)
        # Adding the current and the nonlinear term on the node corresponding to the spine voltage
        U[terminal_node_4+1,t] = U[terminal_node_4+1,t] + Δt*I((t+1)*Δt) + f(U[terminal_node_4+1,t])
        U[:,t+1] .= gmres(A,U[:,t],tol = 1e-8)
        if U[terminal_node_4+1,t+1] < 0
            U[terminal_node_4+1,t+1] = 0
        end
    end
    return U
end

function I(t) 
    if (t < 1) && (t > 0)
        i = 0
    elseif (t >= 1) && (t <= 2)
        i = 10
    elseif (t > 2)
        i = 0
    end
    return i
end

M_t = 20000;
V = zeros(terminal_node_4+2,M_t+1);

V_spine_child_4_tip = solve_linear(A_sparse,V,I);

A = BlockArray(zeros(terminal_node_4+4, terminal_node_4+4), [1, shared_node - 2, 1, terminal_node_1 - shared_node, shared_node_2 - terminal_node_1 - 1, 1, terminal_node_3 - shared_node_2, terminal_node_4 - terminal_node_3,4], [1, shared_node - 2, 1, terminal_node_1 - shared_node, shared_node_2 - terminal_node_1 - 1, 1, terminal_node_3 - shared_node_2, terminal_node_4 - terminal_node_3,4]);

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

A[terminal_node_4+1,terminal_node_4+1] = 0.2*Δt + 1
A[terminal_node_4+1,terminal_node_4+2] = Δt
A[terminal_node_4+2,terminal_node_4+1] = -0.00204*Δt
A[terminal_node_4+2,terminal_node_4+2] = 0.01*Δt + 1;

A[terminal_node_4+3,terminal_node_4+3] = 0.2*Δt + 1
A[terminal_node_4+3,terminal_node_4+4] = Δt
A[terminal_node_4+4,terminal_node_4+3] = -0.00204*Δt
A[terminal_node_4+4,terminal_node_4+4] = 0.01*Δt + 1;

R_spine_mem = 200;

Spine_node = terminal_node_4-1;
Original_value = A[Spine_node,Spine_node]
A[Spine_node,Spine_node] = A[Spine_node,Spine_node] + Δt/(pi*c_m*d5(L)*R_spine_mem);
A[Spine_node,terminal_node_4+1] = -Δt/(pi*c_m*d5(L)*R_spine_mem);

Spine_node_2 = terminal_node_3-1;
Original_value_2 = A[Spine_node_2,Spine_node_2]
A[Spine_node_2,Spine_node_2] = Original_value_2 + Δt/(pi*c_m*d4(L)*R_spine_mem);
A[Spine_node_2,terminal_node_4+3] = -Δt/(pi*c_m*d4(L)*R_spine_mem);

A_sparse = sparse(A);

f(V_Spine_previous) = -Δt/(40^2)*V_Spine_previous^3 + Δt*(0.2+1)/40*V_Spine_previous^2;



function solve_linear(A,U,I)
    for t in 1:(size(U,2)-1)
        # Adding the current and the nonlinear term on the node corresponding to the spine voltage
        U[terminal_node_4+1,t] = U[terminal_node_4+1,t] + Δt*I((t+1)*Δt) + f(U[terminal_node_4+1,t])
        U[terminal_node_4+3,t] = U[terminal_node_4+3,t] + Δt*I((t+1)*Δt) + f(U[terminal_node_4+3,t])
        U[:,t+1] .= gmres(A,U[:,t],tol = 1e-8)
        if U[terminal_node_4+1,t+1] < 0
            U[terminal_node_4+1,t+1] = 0
        end
        if U[terminal_node_4+3,t+1] < 0
            U[terminal_node_4+3,t+1] = 0
        end
    end
    return U
end

function I(t) 
    if (t < 1) && (t > 0)
        i = 0
    elseif (t >= 1) && (t <= 2)
        i = 10
    elseif (t > 2)
        i = 0
    end
    return i
end

M_t = 20000;
V = zeros(terminal_node_4+4,M_t+1);

V_spine_child_4and3_tip = solve_linear(A_sparse,V,I);

t = 0:Δt:(M_t*Δt);

plot(t,V_spine_child_4_tip[100,:],
    xaxis ="time (ms)",yaxis ="voltage (mv)", lw = 2.5,fmt = :svg,label = "individual plateau",linestyle = :dot)
plot!(t,V_spine_child_4and3_tip[100,:],
    xaxis ="time (ms)",yaxis ="voltage (mv)", lw = 2.5,fmt = :svg,label = "two simultaneous plateaus", legendfontsize =8,ylim = (0,2.5))

savefig("FHN and.svg")

d_0 = 0.0001 - 5*0.00001 + 35*0.00001
d_L = 0.00005 - 5*0.00001 + 35*0.00001
d5(x) = ((d_L - d_0)/L)*x + d_0;
d4(x) = ((d_L - d_0)/L)*x + d_0;

A = BlockArray(zeros(terminal_node_4+2, terminal_node_4+2), [1, shared_node - 2, 1, terminal_node_1 - shared_node, shared_node_2 - terminal_node_1 - 1, 1, terminal_node_3 - shared_node_2, terminal_node_4 - terminal_node_3,2], [1, shared_node - 2, 1, terminal_node_1 - shared_node, shared_node_2 - terminal_node_1 - 1, 1, terminal_node_3 - shared_node_2, terminal_node_4 - terminal_node_3,2]);

block_parent = zeros(shared_node-2,shared_node-2);

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

A[terminal_node_4+1,terminal_node_4+1] = 0.2*Δt + 1
A[terminal_node_4+1,terminal_node_4+2] = Δt
A[terminal_node_4+2,terminal_node_4+1] = -0.00204*Δt
A[terminal_node_4+2,terminal_node_4+2] = 0.01*Δt + 1;

R_spine_mem = 200;

Spine_node = terminal_node_4-1;
Original_value = A[Spine_node,Spine_node]
A[Spine_node,Spine_node] = A[Spine_node,Spine_node] + Δt/(pi*c_m*d5(L)*R_spine_mem);
A[Spine_node,terminal_node_4+1] = -Δt/(pi*c_m*d5(L)*R_spine_mem);

A_sparse = sparse(A);

V = zeros(terminal_node_4+2,M_t+1);

function solve_linear(A,U,I)
    for t in 1:(size(U,2)-1)
        # Adding the current and the nonlinear term on the node corresponding to the spine voltage
        U[terminal_node_4+1,t] = U[terminal_node_4+1,t] + Δt*I((t+1)*Δt) + f(U[terminal_node_4+1,t])
        U[:,t+1] .= gmres(A,U[:,t],tol = 1e-8)
        if U[terminal_node_4+1,t+1] < 0
            U[terminal_node_4+1,t+1] = 0
        end
    end
    return U
end

V_spine_child_4_tip = solve_linear(A_sparse,V,I);



A = BlockArray(zeros(terminal_node_4+4, terminal_node_4+4), [1, shared_node - 2, 1, terminal_node_1 - shared_node, shared_node_2 - terminal_node_1 - 1, 1, terminal_node_3 - shared_node_2, terminal_node_4 - terminal_node_3,4], [1, shared_node - 2, 1, terminal_node_1 - shared_node, shared_node_2 - terminal_node_1 - 1, 1, terminal_node_3 - shared_node_2, terminal_node_4 - terminal_node_3,4]);

block_parent = zeros(shared_node-2,shared_node-2);

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

A[terminal_node_4+1,terminal_node_4+1] = 0.2*Δt + 1
A[terminal_node_4+1,terminal_node_4+2] = Δt
A[terminal_node_4+2,terminal_node_4+1] = -0.00204*Δt
A[terminal_node_4+2,terminal_node_4+2] = 0.01*Δt + 1;

A[terminal_node_4+3,terminal_node_4+3] = 0.2*Δt + 1
A[terminal_node_4+3,terminal_node_4+4] = Δt
A[terminal_node_4+4,terminal_node_4+3] = -0.00204*Δt
A[terminal_node_4+4,terminal_node_4+4] = 0.01*Δt + 1;

R_spine_mem = 200;

Spine_node = terminal_node_4-1;
Original_value = A[Spine_node,Spine_node]
A[Spine_node,Spine_node] = A[Spine_node,Spine_node] + Δt/(pi*c_m*d5(L)*R_spine_mem);
A[Spine_node,terminal_node_4+1] = -Δt/(pi*c_m*d5(L)*R_spine_mem);

Spine_node_2 = terminal_node_3-1;
Original_value_2 = A[Spine_node_2,Spine_node_2]
A[Spine_node_2,Spine_node_2] = Original_value_2 + Δt/(pi*c_m*d4(L)*R_spine_mem);
A[Spine_node_2,terminal_node_4+3] = -Δt/(pi*c_m*d4(L)*R_spine_mem);

A_sparse = sparse(A);



function solve_linear(A,U,I)
    for t in 1:(size(U,2)-1)
        # Adding the current and the nonlinear term on the node corresponding to the spine voltage
        U[terminal_node_4+1,t] = U[terminal_node_4+1,t] + Δt*I((t+1)*Δt) + f(U[terminal_node_4+1,t])
        U[terminal_node_4+3,t] = U[terminal_node_4+3,t] + Δt*I((t+1)*Δt) + f(U[terminal_node_4+3,t])
        U[:,t+1] .= gmres(A,U[:,t],tol = 1e-8)
        if U[terminal_node_4+1,t+1] < 0
            U[terminal_node_4+1,t+1] = 0
        end
        if U[terminal_node_4+3,t+1] < 0
            U[terminal_node_4+3,t+1] = 0
        end
    end
    return U
end


V = zeros(terminal_node_4+4,M_t+1);

V_spine_child_4and3_tip = solve_linear(A_sparse,V,I);

t = 0:Δt:(M_t*Δt);
plot(t,V_spine_child_4_tip[100,:],
    xaxis ="time (ms)",yaxis ="voltage (mv)", lw = 2.5,fmt = :svg,label = "individual plateau",linestyle = :dot)
plot!(t,V_spine_child_4and3_tip[100,:],
    xaxis ="time (ms)",yaxis ="voltage (mv)", lw = 2.5,fmt = :svg,label = "two simultaneous plateaus", legendfontsize =8)

savefig("Fig 31 right.svg")

d_avg = zeros(1,35)
amplitudes_indiv = zeros(2,35);

for counter in 1:35
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
    
    d_0 = 0.0001 - 5*0.00001 + counter*0.00001
    d_L = 0.00005 - 5*0.00001 + counter*0.00001
    #d_0_branch_3_4 = 0.0001 #cm - diameter at the begining of the 2st branch
    #d_L_branch_3_4 = 0.00005 #cm - diameter at the end of the 2st branch
    #the linear diameter function for the 2st branch
    d4(x) = ((d_L - d_0)/L)*x + d_0;

    
    d5(x) = ((d_L - d_0)/L)*x + d_0;

    shared_node = Int(round(L/Δx))+1
    terminal_node_1 = shared_node + shared_node - 1
    shared_node_2 = terminal_node_1 + shared_node - 1
    terminal_node_3 = shared_node_2 + shared_node - 1
    terminal_node_4 = terminal_node_3 + shared_node - 1

    A = BlockArray(zeros(terminal_node_4+2, terminal_node_4+2), [1, shared_node - 2, 1, terminal_node_1 - shared_node, shared_node_2 - terminal_node_1 - 1, 1, terminal_node_3 - shared_node_2, terminal_node_4 - terminal_node_3,2], [1, shared_node - 2, 1, terminal_node_1 - shared_node, shared_node_2 - terminal_node_1 - 1, 1, terminal_node_3 - shared_node_2, terminal_node_4 - terminal_node_3,2]);

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
    
    A[terminal_node_4+1,terminal_node_4+1] = 0.2*Δt + 1
    A[terminal_node_4+1,terminal_node_4+2] = Δt
    A[terminal_node_4+2,terminal_node_4+1] = -0.00204*Δt
    A[terminal_node_4+2,terminal_node_4+2] = 0.01*Δt + 1;
    R_spine_mem = 200;
    Spine_node = terminal_node_4-1;
    Original_value = A[Spine_node,Spine_node]
    A[Spine_node,Spine_node] = A[Spine_node,Spine_node] + Δt/(pi*c_m*d5(L)*R_spine_mem);
    A[Spine_node,terminal_node_4+1] = -Δt/(pi*c_m*d5(L)*R_spine_mem);
    A_sparse = sparse(A);
    
    f(V_Spine_previous) = -Δt/(40^2)*V_Spine_previous^3 + Δt*(0.2+1)/40*V_Spine_previous^2;



    function solve_linear(A,U,I)
        for t in 1:(size(U,2)-1)
            # Adding the current and the nonlinear term on the node corresponding to the spine voltage
            U[terminal_node_4+1,t] = U[terminal_node_4+1,t] + Δt*I((t+1)*Δt) + f(U[terminal_node_4+1,t])
            U[:,t+1] .= gmres(A,U[:,t],tol = 1e-8)
            if U[terminal_node_4+1,t+1] < 0
                U[terminal_node_4+1,t+1] = 0
            end
        end
        return U
    end


    function I(t) 
        if (t < 1) && (t > 0)
            i = 0
        elseif (t >= 1) && (t <= 2)
            i = 10
        elseif (t > 2)
            i = 0
        end
        return i
    end

    M_t = 20000;
    V = zeros(terminal_node_4+2,M_t+1);



    V_spine_child_4_tip = solve_linear(A_sparse,V,I);
    
    d_avg[counter] = (d_0 + d_L)/2
    amplitudes_indiv[1,counter] = maximum(V_spine_child_4_tip[1,:])
    amplitudes_indiv[2,counter] = maximum(V_spine_child_4_tip[100,:])
end

d_avg = zeros(1,35)
amplitudes_coop = zeros(2,35);

for counter in 1:35
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

    #d_0_branch_3_4 = 0.0001 #cm - diameter at the begining of the 2st branch
    #d_L_branch_3_4 = 0.00005 #cm - diameter at the end of the 2st branch
    #the linear diameter function for the 2st branch
    #d4(x) = ((d_L_branch_3_4 - d_0_branch_3_4)/L)*x + d_0_branch_3_4;

    d_0 = 0.0001 - 5*0.00001 + counter*0.00001
    d_L = 0.00005 - 5*0.00001 + counter*0.00001
    d4(x) = ((d_L - d_0)/L)*x + d_0;
    d5(x) = ((d_L - d_0)/L)*x + d_0;

    shared_node = Int(round(L/Δx))+1
    terminal_node_1 = shared_node + shared_node - 1
    shared_node_2 = terminal_node_1 + shared_node - 1
    terminal_node_3 = shared_node_2 + shared_node - 1
    terminal_node_4 = terminal_node_3 + shared_node - 1

    A = BlockArray(zeros(terminal_node_4+4, terminal_node_4+4), [1, shared_node - 2, 1, terminal_node_1 - shared_node, shared_node_2 - terminal_node_1 - 1, 1, terminal_node_3 - shared_node_2, terminal_node_4 - terminal_node_3,4], [1, shared_node - 2, 1, terminal_node_1 - shared_node, shared_node_2 - terminal_node_1 - 1, 1, terminal_node_3 - shared_node_2, terminal_node_4 - terminal_node_3,4]);

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
    
    A[terminal_node_4+1,terminal_node_4+1] = 0.2*Δt + 1
    A[terminal_node_4+1,terminal_node_4+2] = Δt
    A[terminal_node_4+2,terminal_node_4+1] = -0.00204*Δt
    A[terminal_node_4+2,terminal_node_4+2] = 0.01*Δt + 1;

    A[terminal_node_4+3,terminal_node_4+3] = 0.2*Δt + 1
    A[terminal_node_4+3,terminal_node_4+4] = Δt
    A[terminal_node_4+4,terminal_node_4+3] = -0.00204*Δt
    A[terminal_node_4+4,terminal_node_4+4] = 0.01*Δt + 1;

    R_spine_mem = 200;

    Spine_node = terminal_node_4-1;
    Original_value = A[Spine_node,Spine_node]
    A[Spine_node,Spine_node] = A[Spine_node,Spine_node] + Δt/(pi*c_m*d5(L)*R_spine_mem);
    A[Spine_node,terminal_node_4+1] = -Δt/(pi*c_m*d5(L)*R_spine_mem);

    Spine_node_2 = terminal_node_3-1;
    Original_value_2 = A[Spine_node_2,Spine_node_2]
    A[Spine_node_2,Spine_node_2] = Original_value_2 + Δt/(pi*c_m*d4(L)*R_spine_mem);
    A[Spine_node_2,terminal_node_4+3] = -Δt/(pi*c_m*d4(L)*R_spine_mem);
    A_sparse = sparse(A);
    
    f(V_Spine_previous) = -Δt/(40^2)*V_Spine_previous^3 + Δt*(0.2+1)/40*V_Spine_previous^2;



    function solve_linear(A,U,I)
        for t in 1:(size(U,2)-1)
            # Adding the current and the nonlinear term on the node corresponding to the spine voltage
            U[terminal_node_4+1,t] = U[terminal_node_4+1,t] + Δt*I((t+1)*Δt) + f(U[terminal_node_4+1,t])
            U[terminal_node_4+3,t] = U[terminal_node_4+3,t] + Δt*I((t+1)*Δt) + f(U[terminal_node_4+3,t])
            U[:,t+1] .= gmres(A,U[:,t],tol = 1e-8)
            if U[terminal_node_4+1,t+1] < 0
                U[terminal_node_4+1,t+1] = 0
            end
            if U[terminal_node_4+3,t+1] < 0
                U[terminal_node_4+3,t+1] = 0
            end
        end
        return U
    end


    function I(t) 
        if (t < 1) && (t > 0)
            i = 0
        elseif (t >= 1) && (t <= 2)
            i = 10
        elseif (t > 2)
            i = 0
        end
        return i
    end

    M_t = 20000;
    V = zeros(terminal_node_4+4,M_t+1);



    V_spine_child_4and3_tip = solve_linear(A_sparse,V,I);
    
    d_avg[counter] = (d_0 + d_L)/2
    amplitudes_coop[1,counter] = maximum(V_spine_child_4and3_tip[1,:])
    amplitudes_coop[2,counter] = maximum(V_spine_child_4and3_tip[100,:])
end

plot(d_avg[:],amplitudes_coop[1,:],seriestype=:scatter,markershape = :square,label = "two simultaneous plateaus (soma)")
plot!(d_avg[:],amplitudes_indiv[1,:],seriestype=:scatter,markershape = :square,linecolor = :red,label = "individual plateau (soma)",legend = :topleft)
plot!(d_avg[:],amplitudes_coop[2,:],seriestype=:scatter,markershape = :star4,label = "two simultaneous plateaus (2nd branch)")
plot!(d_avg[:],amplitudes_indiv[2,:],seriestype=:scatter,markershape = :star4,linecolor = :yellow,label = "individual plateau (2nd branch)",legend = :bottomright,xaxis ="average diameter (cm)",yaxis ="voltage (mv)",xlims = (0,0.0004))

savefig("FHN and diameter.svg")

plot(d_avg[:], amplitudes_coop[1,:]./amplitudes_indiv[1,:],seriestype=:scatter,markershape = :square,label = "2nd branch",xaxis ="average diameter (cm)",xlims = (0,0.0004))
plot!(d_avg[:], amplitudes_coop[2,:]./amplitudes_indiv[2,:],seriestype=:scatter,markershape = :circle,label = "soma",xaxis ="average diameter (cm)",xlims = (0,0.0004),ylim = (1.7,2))

savefig("fig 32 right.svg")
