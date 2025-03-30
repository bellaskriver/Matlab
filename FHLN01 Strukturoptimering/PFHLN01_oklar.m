%Kräver följande .m filer: hooke, plani4e, solveq, extract, designDomain, plani4s, eldraw2
XVAL = projektFHLN01;

%Huvudfunktion
function XVAL = projektFHLN01
format short; clc; close all;

%Inparametrar mesh
Lx = 0.06;               
Ly = 0.02;               
le = 0.0005;              
nelm = (Lx/le) * (Ly/le);
[coord, ~, enod, edof, Ex, Ey] = designDomain(Lx, Ly, le);
    
%Inparametrar upplag
bc_xmin = 0.0085;  
bc_xmax = 0.0105;  
%bcupp = [bc_xmin, bc_xmax]; BELLA

%Inparametrar spänningsvillkor
p_stress = 2;
q_stress = 3;
sigma_allow = 1e6; %Maxiaml tillåten spänning (Pa)

%Inparametrar
t = 0.01; %Tjocklek (m)
de = 1e-9; %Resudual styvhet delta
kraft = -50; %Kraft (N)  
xval = 0.25 * ones(nelm, 1); %Designvariabler
r_min_compli = 3 * le;  
xforce_left = 2 * le;

%Inparametrar Heaviside
beta = 1;
eta = 0.5;
    
%Skapa filtermatris
M = filterDesignVariable(r_min_compli, coord, enod, nelm);
    
%Inparametrar MMA optimering
m = 2; %m=1 volymsvillkor, m=2 volymsvillkor + spänningsvillkor
maxIter = 70; %Antal iterationer
xmin = zeros(nelm, 1); 
xmax = ones(nelm, 1);
xold1 = xval;    
xold2 = xval;
low = xmin;      
upp = xmax;
a0 = 1;
a = zeros(m,1);
c = nelm * ones(m,1);
d = ones(m,1);
f0_history = zeros(maxIter,1); 

%Optimeringsloop
for iter = 1:maxIter
    fprintf('\nIteration %d/%d\n', iter, maxIter);

    %Öka SIMP penaliseringen linjärt
    p = 3 + (iter-1)/(maxIter-1)*(5-3); 
        
    xval_filtered = M * xval; %Filtrerade designvariabeler
    xval_projected = projectDesign(xval_filtered, beta, eta); %Heaviside projekterade designvariabeler
    dproj = derivativeProject(xval_filtered, beta, eta); %Derivatan av Heaviside projektion
        
    %Lös FEM
    [kelem, K, ac, f, Ed, bc, D] = FEMlosare(edof, Ex, Ey, xval_projected, nelm, Lx, Ly, le, de, p, xforce_left, kraft, bc_xmin, bc_xmax);
        
    %Målfunktion: komplians
    [dc_vec, C] = adjoint_compli(kelem, f, Ed, nelm, xval_projected, ac, p, de, M);
    dc_vec = dc_vec .* dproj;
    f0val = C; %Komplians
    df0dx = dc_vec; %Derivatan av komplians
        
    %Volymsvillkor
    volFrac = sum(xval_projected)/nelm;
    g_vol = volFrac/0.25 - 1; %Volymsvillkor
    dfdx_vol = (1/nelm) * dproj'; %Derivatan av volymsvillkor
        
    %Spänningsvillkor
    if m == 2
        [g_stress, dg_stress] = adjoint_stress(edof, Ex, Ey, kelem, K, ac, bc, D, xval_projected, p, sigma_allow, de, M, nelm);
        dg_stress = (dg_stress .* dproj);
        fval = [g_vol; g_stress];
        dfdx = [dfdx_vol; dg_stress(:)'];
    else
        fval = g_vol;
        dfdx = dfdx_vol;
    end
        
    %Kör mmasub
    [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low, upp] = mmasub(m, nelm, iter, xval, xmin, xmax, xold1, xold2, f0val, df0dx, fval, dfdx, low, upp, a0, a, c, d);
   
    %Spara resultat av iteration
    xval = xmma;
    f0_history(iter) = f0val;
    xold2 = xold1;
    xold1 = xval;
        
    %Visa resultat av iteration
    figure(1); clf;
    patch(Ex', Ey', xval_projected', 'EdgeColor','none');
    colorbar; axis equal tight;
    set(colorbar, 'Limits',[0 1] );
    title(['Iteration ', num2str(iter)]);
    drawnow;
        
    %Öka beta var 10:e iteration
    if mod(iter,10)==0
        beta = beta * 1.5;
    end
end

%Slutligt resultat
XVAL = xval_projected;

%Visa kraftplacering
[row_indices, ~, values] = find(f);
fprintf('\nFrihetsgrad\tKraft\n');
fprintf('%d\t\t%.6f\n', [row_indices, values]');
    
% %Plotta resultat LUCAS
% k = 1; 
% kt = 1; 
% kk = 1; 
% kkk = 1; 
% kkkk = 1; 
% kkkkk = 1;
% nx = round(Lx/le);
% ny = round(Ly/le);
% xmat = reshape(XVAL, [nx, ny]);
% meshplot(Ex, Ey, Ed, k, kt, kk, kkk, kkkk, kkkkk, XVAL, coord, dof, Lx, le, [], maxIter, f0_history, xmat);

%Känslighetsanalys
elementToPerturb = 89; %Element som undersöks
compareStressSens(elementToPerturb, xval, Lx, Ly, le, de, p_stress, sigma_allow, q_stress, bc_xmin, bc_xmax, xforce_left, kraft, M, beta, eta, edof, Ex, Ey, kelem, ac, bc, D, nelm);

%Skapa STL fil
generate_stl(0.5, 0.01, 'myDesign.stl', XVAL, coord, enod, Lx);
fprintf('\nSTL file "myDesign.stl" was created.\n');
end

%Underfunktioner
%FEM lösare
function [kelem, K, ac, f, Ed, bc, D] = FEMlosare(edof, Ex, Ey, xval_tilde, nelm, Lx, Ly, le, de, p, xforce_left, kraft, bc_xmin, bc_xmax)
ptype = 1; E = 4.5e9; nu = 0.32; 
D = hooke(ptype, E, nu);
ir = 2;
ep = [ones(nelm,1), xval_tilde, ir*ones(nelm,1)];

maxedof = max(edof(:));
K = sparse(maxedof, maxedof);
f = sparse(maxedof, 1);

bc = platt(Lx, Ly, le, bc_xmin, bc_xmax);

force_x_min = Lx - xforce_left;
force_x_max = Lx;
i_x_min_force = round(force_x_min/le)+1;
i_x_max_force = round(force_x_max/le)+1;

num_nodes_x = Lx/le + 1;
num_nodes_y = Ly/le + 1;
top_edge_start = (num_nodes_y - 1)*num_nodes_x + 1;
node_indices_x_range = (top_edge_start + i_x_min_force - 1 : top_edge_start + i_x_max_force - 1)';
DOF_force_y = 2*node_indices_x_range;
f(DOF_force_y) = kraft / length(DOF_force_y);

kelem = cell(nelm,1);
for i = 1:nelm
    Ke = plani4e(Ex(i,:), Ey(i,:), ep(i,:), D);
    kelem{i} = Ke;
    idx = edof(i,2:end);
    penalVal = de + (1-de) * xval_tilde(i)^p;
    K(idx,idx) = K(idx,idx) + Ke * penalVal; %Styvhetsmatris
    end
[ac, ~] = solveq(K, f, bc); %Global förskjutningsvektor
Ed = extract(edof, ac); %Elementförskjutingar
end

%Upplagsvillkor
function [bc] = platt(Lx, Ly, le, bc_xmin, bc_xmax)
    i_x_min = round(bc_xmin/le) + 1;
    i_x_max = round(bc_xmax/le) + 1;
    node_indices_x_range = i_x_min:i_x_max;
    DOF_x = 2*node_indices_x_range - 1;
    DOF_y = 2*node_indices_x_range;
    bcupp = [sort(DOF_y(:)), zeros(length(DOF_y),1)];
    numPtsY = Ly/le + 1;
    bchvec = [];
    for i = 1:numPtsY
        bch = ((Lx/le)*2 + 1)*i + i - 1;
        bchvec = [bchvec, bch];
    end
    bc_right_edge = [bchvec(:), zeros(length(bchvec),1)];
    bc = sortrows([bc_right_edge; bcupp]);
end

%Complians med adjoint metod
function [dc_vec, C] = adjoint_compli(kelem, f, Ed, nelm, xval_projected, ac, p, de, M)
C = f' * ac;
dc_vec = zeros(nelm, 1);
for i = 1:nelm
    fac = de + (1 - de) * xval_projected(i)^p;
    dc_vec(i) = -p * (1-de) * xval_projected(i)^(p-1) * Ed(i,:) * kelem{i} * Ed(i,:)';
end

%Normalisera
C = C / 40;
dc_vec = (dc_vec / 40);

%Filtrera
dc_vec = M * dc_vec;
end

%Spänningsvillkor med adjoint metod
function [g_stress, dg_stress] = adjoint_stress(edof, Ex, Ey, kelem, K, ac, bc, D, xval_projected, p, sigma_allow, de, M, nelm)
ndof = size(K,1);
vm = zeros(nelm,1);
P = [2 -1 0; -1 2 0; 0 0 6];
for e = 1:nelm
    B = calcB4nodeCenter(Ex(e,:), Ey(e,:));
    edof_e = edof(e,2:end);
    u_e = ac(edof_e);
    sigma_e = D * B * u_e;
    vm(e) = sqrt(0.5 * (sigma_e' * P * sigma_e));
end
vm_p = vm.^p;
Phi = ((1/nelm) * sum(vm_p))^(1/p);
g_stress = (Phi / sigma_allow) - 1;
R = zeros(ndof,1);
preFactor = (1/p) * Phi^(1-p) / nelm;
for e = 1:nelm
    if vm(e) < 1e-14, continue; end
        fac_e = p * vm(e)^(p-1);
        B = calcB4nodeCenter(Ex(e,:), Ey(e,:));
        edof_e = edof(e,2:end);
        u_e = ac(edof_e);
        sigma_e = D * B * u_e;
        grad_vm = (sigma_e' * P / vm(e)) * (D * B);
        R(edof_e) = R(edof_e) + fac_e * grad_vm';
    end
    R = preFactor * R;
    free_dofs = setdiff(1:ndof, bc(:,1));
    lambda_stress = zeros(ndof,1);
    K_reduced = K(free_dofs, free_dofs);
    lambda_stress(free_dofs) = K_reduced \ (-R(free_dofs));
    U = ac;
    dgPhi_dx = zeros(nelm,1);
    for e = 1:nelm
        edof_e = edof(e,2:end);
        Ue = U(edof_e);
        penal_deriv = p * (1-de) * xval_projected(e)^(p-1);
        dK_dx_e = penal_deriv * kelem{e};
        dgPhi_dx(e) = lambda_stress(edof_e)' * (dK_dx_e * Ue);
    end
dg_stress = (1 / sigma_allow) * dgPhi_dx;

%Normalisera
g_stress = g_stress / 40;
dg_stress = dg_stress / 40;

%Filtrera
dg_stress = M * dg_stress;
end

%B matris för spänningsvillkor
function B = calcB4nodeCenter(ex, ey)
xi = 0.0;
eta = 0.0;
dNdxi = 0.25 * [-(1-eta), (1-eta), (1+eta), -(1+eta); -(1-xi), -(1+xi), (1+xi), (1-xi)];
x = ex(:);
y = ey(:);
J = dNdxi * [x, y];
invJ = inv(J);
dN = invJ * dNdxi;
B = zeros(3,8);
for i = 1:4
    dNdx = dN(1,i);
    dNdy = dN(2,i);
    B(1,2*i-1) = dNdx;
    B(2,2*i)   = dNdy;
    B(3,2*i-1) = dNdy;
    B(3,2*i)   = dNdx;
end
end

%Filtermatris
function M = filterDesignVariable(r_o, coord, enod, nelm)
    centroids = zeros(nelm, 2);
    for e = 1:nelm
        node_indices = unique(enod(e,2:end));
        node_indices = node_indices(node_indices>0 & node_indices<=size(coord,1));
        centroids(e,:) = mean(coord(node_indices,:),1);
    end
    M = zeros(nelm, nelm);
    for e = 1:nelm
        c_e = centroids(e,:);
        sumW = 0;
        for j = 1:nelm
            dist_ej = norm(c_e - centroids(j,:));
            if dist_ej < r_o
                sumW = sumW + (r_o - dist_ej);
            end
        end
        for j = 1:nelm
            dist_ej = norm(c_e - centroids(j,:));
            if dist_ej < r_o
                M(e,j) = (r_o - dist_ej) / max(sumW, 1e-15);
            end
        end
    end
end

%Heavyside projektion funktion
function xproj = projectDesign(rho_f, beta, eta)
    top = tanh(beta*eta) + tanh(beta*(rho_f - eta));
    bot = tanh(beta*eta) + tanh(beta*(1 - eta));
    xproj = top ./ bot;
end

%Heavyside projektion derivata
function dproj = derivativeProject(rho_f, beta, eta)
    bot = tanh(beta*eta) + tanh(beta*(1 - eta));
    dproj = (beta * sech(beta*(rho_f - eta)).^2) ./ bot;
end

%Känsligheter
function [g_stress, dg_stress] = evaluateStressConstraint(xval_in, Lx, Ly, le, de, p_stress, sigma_allow, M, beta, eta)
[coord, dof, enod, edof, Ex, Ey] = designDomain(Lx, Ly, le);
nelm = (Lx/le)*(Ly/le);
xval_filtered = M * xval_in;
xval_projected = projectDesign(xval_filtered, beta, eta);
dproj = derivativeProject(xval_filtered, beta, eta);
    
xforce_left_fixed = 2 * le;
kraft_fixed = -50;
bc_xmin_fixed = 0.0085;
bc_xmax_fixed = 0.0105;
[kelem, K, ac, f, Ed, bc, D] = FEMlosare(edof, Ex, Ey, xval_projected, nelm, Lx, Ly, le, de, p_stress, xforce_left_fixed, kraft_fixed, bc_xmin_fixed, bc_xmax_fixed);
    
[g_stress, dg_stress] = adjoint_stress(edof, Ex, Ey, kelem, K, ac, bc, D, xval_projected, p_stress, sigma_allow, de, M, nelm);  
dg_stress = dg_stress .* dproj;
end

%Jämför känsligheter med central difference method
function compareStressSens(elementToPerturb, xval_baseline, Lx, Ly, le, de, p_stress, sigma_allow, q_stress, bc_xmin, bc_xmax, xforce_left, kraft, M, beta, eta, edof, Ex, Ey, kelem, ac, bc, D, nelm)
h = 1e-6; % Steglängd för perturbation
xval_perturbed = xval_baseline;
xval_perturbed(elementToPerturb) = xval_perturbed(elementToPerturb) + h;

%Numerisk känslighet, baseline vs pertuberad
[g_stress_base, dg_stress_base] = evaluateStressConstraint(xval_baseline, Lx, Ly, le, de, p_stress, sigma_allow, M, beta, eta);
[g_stress_pert, ~] = evaluateStressConstraint(xval_perturbed, Lx, Ly, le, de, p_stress, sigma_allow, M, beta, eta);
numSens = (g_stress_pert - g_stress_base) / h;
    
%Analytisk känslighet
analyticSens = dg_stress_base(elementToPerturb);

%Jämförelse
ratio = numSens / analyticSens;
fprintf('\n--- Sensitivity Check for Element %d ---\n', elementToPerturb);
fprintf('Numerical  = %12.6e\n', numSens);
fprintf('Analytical = %12.6e\n', analyticSens);
fprintf('Ratio      = %12.6e\n', ratio);
end

%function meshplot(Ex, Ey, Ed, k, kt, kk, kkk, kkkk, kkkkk, XVAL, coord, dof, Lx, le, vonMises_effective, maxIter, f0_history, xmat) LUCAS
%     if k==1
%         figure();
%         patch(Ex', Ey', 1);
%         plotpar=[1 2 0];
%         eldraw2(Ex, Ey, plotpar);
%         if kt==1
%             text(coord(:,1), coord(:,2), num2str(dof(:,1)))
%         end
%         axis equal;
%     end
%     if kk==1
%         figure();
%         hold on; axis equal;
%         sfac = 500;
%         for i = 1:length(Ex)
%             x_orig = Ex(i,:);
%             y_orig = Ey(i,:);
%             ux = Ed(i,[1,3,5,7]);
%             uy = Ed(i,[2,4,6,8]);
%             x_def = x_orig + sfac*ux;
%             y_def = y_orig + sfac*uy;
%             plot([x_orig, x_orig(1)], [y_orig, y_orig(1)], 'g-', 'LineWidth', 0.5);
%             plot([x_def, x_def(1)], [y_def, y_def(1)], 'r-', 'LineWidth', 1);
%         end
%         title('Odeformerad vs. deformerad');
%         hold off;
%     end
%     if kkk==1
%         figure();
%         xmat_flipped = flip(xmat);
%         imagesc(xmat_flipped);
%         colormap(flipud(hot));
%         colorbar;
%         clim([0 1]);
%         axis equal;
%         title('Densitetsfördelning');
%     end  
%     if kkkk==1
%         figure();
%         patch(Ex', Ey', vonMises_effective');
%         colormap('jet');
%         colorbar;
%         title('Von Mises Stress Distribution');
%         axis equal;
%     end
%     if kkkkk==1
%         figure();
%         plot(1:maxIter, f0_history, 'b-o', 'LineWidth', 1.5);
%         xlabel('Iteration');
%         ylabel('Complians');
%         title('Convergens av complians');
%         grid on;
%     end
% end

%mmasub (modifierad: move=0.05)
function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
%
%    This function mmasub performs one MMA-iteration, aimed at
%    solving the nonlinear programming problem:
%         
%      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
%    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m
%                xmin_j <= x_j <= xmax_j,    j = 1,...,n
%                z >= 0,   y_i >= 0,         i = 1,...,m
%*** INPUT:
%
%   m    = The number of general constraints.
%   n    = The number of variables x_j.
%  iter  = Current iteration number ( =1 the first time mmasub is called).
%  xval  = Column vector with the current values of the variables x_j.
%  xmin  = Column vector with the lower bounds for the variables x_j.
%  xmax  = Column vector with the upper bounds for the variables x_j.
%  xold1 = xval, one iteration ago (provided that iter>1).
%  xold2 = xval, two iterations ago (provided that iter>2).
%  f0val = The value of the objective function f_0 at xval.
%  df0dx = Column vector with the derivatives of the objective function
%          f_0 with respect to the variables x_j, calculated at xval.
%  fval  = Column vector with the values of the constraint functions f_i,
%          calculated at xval.
%  dfdx  = (m x n)-matrix with the derivatives of the constraint functions
%          f_i with respect to the variables x_j, calculated at xval.
%          dfdx(i,j) = the derivative of f_i with respect to x_j.
%  low   = Column vector with the lower asymptotes from the previous
%          iteration (provided that iter>1).
%  upp   = Column vector with the upper asymptotes from the previous
%          iteration (provided that iter>1).
%  a0    = The constants a_0 in the term a_0*z.
%  a     = Column vector with the constants a_i in the terms a_i*z.
%  c     = Column vector with the constants c_i in the terms c_i*y_i.
%  d     = Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
%     
%*** OUTPUT:
%
%  xmma  = Column vector with the optimal values of the variables x_j
%          in the current MMA subproblem.
%  ymma  = Column vector with the optimal values of the variables y_i
%          in the current MMA subproblem.
%  zmma  = Scalar with the optimal value of the variable z
%          in the current MMA subproblem.
%  lam   = Lagrange multipliers for the m general MMA constraints.
%  xsi   = Lagrange multipliers for the n constraints alfa_j - x_j <= 0.
%  eta   = Lagrange multipliers for the n constraints x_j - beta_j <= 0.
%   mu   = Lagrange multipliers for the m constraints -y_i <= 0.
%  zet   = Lagrange multiplier for the single constraint -z <= 0.
%   s    = Slack variables for the m general MMA constraints.
%  low   = Column vector with the lower asymptotes, calculated and used
%          in the current MMA subproblem.
%  upp   = Column vector with the upper asymptotes, calculated and used
%          in the current MMA subproblem.
%
%epsimin = sqrt(m+n)*10^(-9);
epsimin = 10^(-7);
raa0 = 0.00001;
move = 0.05; %Modifierad för att underlätta konvergens
albefa = 0.1;
asyinit = 0.5;
asyincr = 1.2;
asydecr = 0.7;
eeen = ones(n,1);
eeem = ones(m,1);
zeron = zeros(n,1);

% Calculation of the asymptotes low and upp :
if iter < 2.5
  low = xval - asyinit*(xmax-xmin);
  upp = xval + asyinit*(xmax-xmin);
else
  zzz = (xval-xold1).*(xold1-xold2);
  factor = eeen;
  factor(find(zzz > 0)) = asyincr;
  factor(find(zzz < 0)) = asydecr;
  low = xval - factor.*(xold1 - low);
  upp = xval + factor.*(upp - xold1);
  lowmin = xval - 10*(xmax-xmin);
  lowmax = xval - 0.01*(xmax-xmin);
  uppmin = xval + 0.01*(xmax-xmin);
  uppmax = xval + 10*(xmax-xmin);
  low = max(low,lowmin);
  low = min(low,lowmax);
  upp = min(upp,uppmax);
  upp = max(upp,uppmin);
end

% Calculation of the bounds alfa and beta :

zzz1 = low + albefa*(xval-low);
zzz2 = xval - move*(xmax-xmin);
zzz  = max(zzz1,zzz2);
alfa = max(zzz,xmin);
zzz1 = upp - albefa*(upp-xval);
zzz2 = xval + move*(xmax-xmin);
zzz  = min(zzz1,zzz2);
beta = min(zzz,xmax);

% Calculations of p0, q0, P, Q and b.

xmami = xmax-xmin;
xmamieps = 0.00001*eeen;
xmami = max(xmami,xmamieps);
xmamiinv = eeen./xmami;
ux1 = upp-xval;
ux2 = ux1.*ux1;
xl1 = xval-low;
xl2 = xl1.*xl1;
uxinv = eeen./ux1;
xlinv = eeen./xl1;
%
p0 = zeron;
q0 = zeron;
p0 = max(df0dx,0);
q0 = max(-df0dx,0);
pq0 = 0.001*(p0 + q0) + raa0*xmamiinv;
p0 = p0 + pq0;
q0 = q0 + pq0;
p0 = p0.*ux2;
q0 = q0.*xl2;
%
P = sparse(m,n);
Q = sparse(m,n);
P = max(dfdx,0);
Q = max(-dfdx,0);
PQ = 0.001*(P + Q) + raa0*eeem*xmamiinv';
P = P + PQ;
Q = Q + PQ;
P = P * spdiags(ux2,0,n,n);
Q = Q * spdiags(xl2,0,n,n);
b = P*uxinv + Q*xlinv - fval ;
%
%%% Solving the subproblem by a primal-dual Newton method
[xmma,ymma,zmma,lam,xsi,eta,mu,zet,s] = ...
subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);
end

%generate_stl (modifierad: symmetri kring höger sida)
function generate_stl(threshold_den_val, thickness, name, xval_tilde, coord, enod, Lx)
el = [];
for i = 1:length(xval_tilde)
    if xval_tilde(i) >= threshold_den_val
        el = [el, i];
    end
end
nodes = [];
for i = el
    nodes = [nodes; enod(i,2:5)'];
end
node_num = unique(nodes);
X = coord(node_num,1);
Y = coord(node_num,2);
Z = zeros(length(X),1);
Z = [Z; thickness.*ones(length(X), 1)];
X = [X; X];
Y = [Y; Y];
X = X - Lx; %Vänster halva har betraktats
X = [X; -X];
Y = [Y; Y];
Z = [Z; Z];
shp = alphaShape(X, Y, Z);
[bf,P] = boundaryFacets(shp);
stlwrite(triangulation(bf, P), name);
end

