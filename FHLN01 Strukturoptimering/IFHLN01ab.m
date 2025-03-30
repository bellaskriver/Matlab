%Kräver följande .m filer: hooke, plani4e, solveq, extract, designDomain, plani4s, eldraw2
projektFHLN01;

%Huvudfunktion
function projektFHLN01
format short; clc; close all;

%Inparametrar mesh
Lx = 0.06;
Ly = 0.02;
le = 0.001;
nelm = (Lx/le) * (Ly/le);
[coord, dof, enod, edof, Ex, Ey] = designDomain(Lx, Ly, le);

%Inparametrar upplag
bc_xmin = 0.0085;
bc_xmax = 0.0105;

%Inparametrar
t = 0.15; %Tjocklek (m)
de = 1e-9; %Resudual styvhet delta
kraft = -50; %Kraft (N)
xval = t * ones(nelm, 1); %Designvariabler
r_min_compli = 4 * le; %Filterradie
xforce_left = 2 * le;

%Skapa filtermatris
M = filterDesignVariable(r_min_compli, coord, enod, nelm);

%Inparametrar MMA optimering
m = 1; %m=1 volymsvillkor, m=2 volymsvillkor + spänningsvillkor
maxIter = 150; %Antal iterationer
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

xval_filtered = M * xval;

%Öka SIMP penaliseringen linjärt
p = 3 + (iter-1)/(maxIter-1)*(5-3);

%Lös FEM
[kelem, K, ac, f, Ed, bc, D] = FEMlosare(edof, Ex, Ey, xval, nelm, Lx, Ly, le, de, p, xforce_left, kraft, bc_xmin, bc_xmax);

%Målfunktion: komplians
[dc_vec, C] = adjoint_compli(kelem, f, Ed, nelm, xval, ac, p, de);
currVol=sum(xval)*(le^2);
maxVol=sum(nelm)*(le^2);

%Filtrera
dc_vec = M * dc_vec;

f0val=C;
df0dx=dc_vec;

volFrac = sum(M*xval)/nelm;
fval = volFrac/0.25-1;   % (1 st constraint)
dfdx=ones(m,nelm)*M*(4*(1/nelm));

%Kör mmasub
[xmma, ~, ~, ~, ~, ~, ~, ~, ~, low, upp] = mmasub(m, nelm, iter, xval, xmin, xmax, xold1, xold2, f0val, df0dx, fval, dfdx, low, upp, a0, a, c, d);

%Spara resultat av iteration
xval = xmma;
f0_history(iter) = f0val;
xold2 = xold1;
xold1 = xval;
figure(1); clf;
patch(Ex', Ey', xval', 'EdgeColor', 'none');
colorbar;

% Set the axis to be equal and tight
axis equal tight;

% Create a gray colormap where 0 is white and 1 is black, with 10 levels
colormap(repmat(linspace(1, 0, 10)', 1, 3)); % 10-level grayscale (0: white, 1: black)

% Set color axis limits
caxis([0 1]);

% Manually set the colorbar ticks
colorbar;
set(colorbar, 'Ticks', 0:0.1:1); % Ticks at intervals of 0.1

% Add title with iteration number
title(['Iteration ', num2str(iter), ' : Filterradie = ', num2str(r_min_compli/le), ' element']);

% Update the plot
drawnow;
end

%Slutligt resultat
XVAL = xval;
disp(C)
save("xfile","xval")

%Visa kraftplacering
[row_indices, ~, values] = find(f);
fprintf('\nFrihetsgrad\tKraft\n');
fprintf('%d\t\t%.6f\n', [row_indices, values]');

generate_stl(0.5, 0.01, 'myDesign.stl', XVAL, coord, enod, Lx);
fprintf('\nSTL file "myDesign.stl" was created.\n');

%Plotta resultat
k = 0;
kt = 0;
kk = 0;
kkk = 1;
kkkk = 1;
kkkkk = 1;
nx = round(Lx/le);
ny = round(Ly/le);
xmat = reshape(XVAL, [nx, ny]);
meshplot(edof, Ex, Ey, Ed, k, kt, kk, kkk, kkkk, kkkkk, coord, dof, maxIter, f0_history, xmat)
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
function [dc_vec, C] = adjoint_compli(kelem, f, Ed, nelm, xval, ac, p, de);
C = f' * ac;
dc_vec = zeros(nelm, 1);

for i = 1:nelm
dc_vec(i) = -p * (1-de) * xval(i)^(p-1) * Ed(i,:) * kelem{i} * Ed(i,:)';
end

%Normalisera
C=C*1000;
dc_vec = (dc_vec*1000);

end

%Plottar resultat
function meshplot(edof, Ex, Ey, Ed, k, kt, kk, kkk, kkkk, kkkkk, coord, dof, maxIter, f0_history, xmat)
if k==1
figure();
patch(Ex', Ey', 1);
plotpar=[1 2 0];
eldraw2(Ex, Ey, plotpar);
if kt==1
text(coord(:,1), coord(:,2), num2str(dof(:,1)))
end
axis equal;
end

if kk==1
figure();
hold on; axis equal;
sfac = 10;
for i = 1:length(Ex)
x_orig = Ex(i,:);
y_orig = Ey(i,:);
ux = Ed(i,[1,3,5,7]);
uy = Ed(i,[2,4,6,8]);
x_def = x_orig + sfac*ux;
y_def = y_orig + sfac*uy;
plot([x_orig, x_orig(1)], [y_orig, y_orig(1)], 'g-', 'LineWidth', 0.5);
plot([x_def, x_def(1)], [y_def, y_def(1)], 'r-', 'LineWidth', 1);
end
title('Odeformerad vs. deformerad');
hold off;
end

if kkk==1
figure()
xmat_rotated = rot90(xmat, -3); 
imagesc(xmat_rotated);
colormap(flipud(gray)); 
colorbar;
clim([0 1]);
axis equal;
title('Densitetsfördelning');
end

if kkkkk==1
figure();
plot(1:maxIter, f0_history, 'b-o', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Komplians');
title('Konvergens av komplians');
grid on;
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
move = 0.05;
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

%Mesh
function [coord, dof, enod, edof, Ex, Ey] = designDomain(Lx, Ly, le); %Skapar mesh
% Generates a recangular design domain with length Lx in the x-direction
% and Ly in the y-direction. Four node elements are used with the side
% length le (same in both x- and y-direction). You can check the mesh by
% running patch(Ex', Ey', 1) in the command window.
% Note that symmetry is considered along x=0. The symmetry boundary conditions
% are imposed via the matrix bc which contains the degrees of freedom where
% displacements are prescribed in the first column, and the value in the
% seconed column. Extend this matrix with other desirable boundary conditions.
elem_x = Lx/le; elem_y = Ly/le; nelm = elem_x*elem_y;
nodes_x = elem_x + 1; nodes_y = elem_y + 1; nnod = nodes_x*nodes_y;
%coord
coord = zeros(nnod, 2);
node = 1;
for y = 0:nodes_y-1
    for x = 0:nodes_x-1
        coord(node,:) = [x*le y*le];
        node = node + 1;
    end
end
%coord done
%dof
dof = zeros(nnod, 2);
for i = 1:nnod
    dof(i,:) = [i*2-1 i*2];
end
%dof done
%enod
enod = zeros(nelm, 5);
enod(:,1) = 1:nelm;
enod(1,2:5) = [1 2 nodes_x+2 nodes_x+1];
for i = 2:nelm
    if (mod(i, elem_x) == 1)
        enod(i, 2:5) = enod(i-1, 2:5) + 2;
    else
        enod(i, 2:5) = enod(i-1, 2:5) + 1;
    end
end
%enod done
%Ex, Ey
Ex = zeros(nelm, 4);
Ey = zeros(nelm, 4);
for i = 1:nelm
    Ex(i,:) = coord(enod(i, 2:5), 1);
    Ey(i,:) = coord(enod(i, 2:5), 2);
end
%Ex, Ey done
%edof
edof = zeros(nelm, 9);
edof(:,1) = 1:nelm;
for i = 1:nelm
    edof(i, 2:9) = [enod(i, 2)*2-1 enod(i, 2)*2 enod(i, 3)*2-1 enod(i, 3)*2 enod(i,4)*2-1 enod(i, 4)*2 enod(i, 5)*2-1 enod(i, 5)*2];
end
%edof done
%symmetry bc
bc = [];
for i = 1:nnod
    if (coord(i, 1) == 0.0)
        bc = [bc ; i*2 - 1 0];
    end
end
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
