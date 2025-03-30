%Assignment 2 FEM Flöde - Isabella McAuley Skriver & Andrea Nilsson

%% Utan lerlager

clear all
close all
format compact
format short

%Punkter
Vertices = [
    0 0 %1
    10 0 %2
    20 0 %3
    30 0 %4
    40 0 %5
    0 3 %6
    10 3 %7
    20 3 %8
    30 3 %9
    40 3 %10
    0 6 %11
    10 6 %12
    20 6 %13
    30 6 %14
    40 6 %15
    0 10 %16
    10 10 %17
    20 10 %18
    30-0.025 10 %19
    30+0.025 10 %20
    40 10 %21
    ];

%Linjer
Segments = [
    1 2 %1
    2 3 %2
    3 4 %3
    4 5 %4
    6 7 %5
    7 8 %6
    8 9 %7
    9 10 %8
    11 12 %9
    12 13 %10
    13 14 %11
    14 15 %12
    16 17 %13
    17 18 %14
    18 19 %15
    20 21 %16
    1 6 %17
    6 11 %18
    11 16 %19
    2 7 %20
    7 12 %21
    12 17 %22
    3 8 %23
    8 13 %24
    13 18 %25
    4 9 %36
    9 14 %27
    14 19 %28
    14 20 %29
    5 10 %30
    10 15 %31
    15 21 %32
    ];

%Ytor
Surfaces=[
    1 20 5 17 %1
    2 23 6 20 %2
    3 26 7 23 %3
    4 30 8 26 %4
    5 21 9 18 %5
    6 24 10 21 %6
    7 27 11 24 %7
    8 31 12 27 %8
    9 22 13 19 %9
    10 25 14 22 %10
    11 28 15 25 %11
    12 32 16 29 %12
    ];

%Antal indelningar per linje
Seed =[10 10 20 20 10 10 20 20 10 10 20 20 10 10 20 20 3 6 8 3 6 8 3 6 8 3 6 8 8 3 6 8];

segp=[Seed];
nen=4;
dofsPerNode=1;
mp=[dofsPerNode, nen];

%Illustrerar elementindelning
geomdraw2(Vertices,Segments,Surfaces,segp,mp)

%Generera mesh
[Coord Edof Dof meshdb ]=strMeshgen(Vertices,Segments,Surfaces,segp,mp);

%Generera elementkoordinater
[Ex,Ey]=coordxtr(Edof,Coord,Dof,nen);

%Rita mesh med numrerade element
figure(2)
eldraw2(Ex,Ey,[1 4 1],Edof(:,1));

%Rita mesh med numrerade frihetsgrader
figure(3)
eldraw2(Ex,Ey,[1 4 0]);
text(Coord(:,1),Coord(:,2),num2str(Dof(:)))

%Skapa elementstyvhet, assemblera global styvhetsmatris, bestäm randvilkor
%och lös ekvationssystemet

%Egenskaper
ep=[1];
D=eye(2)*10^-4;
nel=size(Edof,1);
ndof=max(Edof(:));

%Nollor tas bort, elementstyvhetsmatris och assemblering
K=sparse(ndof,ndof);
f=sparse(ndof,1);
for i=1:nel
    Ke=flw2qe(Ex(i,:),Ey(i,:),ep,D);
    K=sparse_assem(Edof(i,:),K,Ke);
end

%Randvillkor
bc1=extrSeg([13 14 15]',meshdb,1);
bc2=extrSeg([16]',meshdb,1);

%Definera randvillkor
bc=[bc1, ones(length(bc1),1)*7; bc2, ones(length(bc2),1)*0];

%Lös ekvationssystemet
[a,r]=solveq(K,f,bc);
Ed=extract_ed(Edof,a);

%Beräkna element flux och gradient
for i=1:nel
    [Es(i,:),Et(i,:)]=flw2qs(Ex(i,:),Ey(i,:),ep,D,Ed(i,:));
end

%Plotta flöde
figure(4)
eldraw2(Ex,Ey,[1 4 0]);
[sfac]=elflux2(Ex(350,:),Ey(350,:),Es(350,:));
elflux2(Ex,Ey,Es,[1 2],sfac);

%Plotta tryck
figure(5)
eldraw2(Ex,Ey, [1 4 0]);
eliso2(Ex,Ey,Ed,10,[1 2 1]);

%Plotta färger
figure(6)
fill(Ex',Ey',Ed')
axis('equal')

%Flöde = 3.9566e-04
flodein=sum(r(bc1))
flodeut=sum(r(bc2))

%Gradient risk för erosion
[Gradient]=((Et(:,2).^2+Et(:,1).^2).^0.5);
Gradkrit=Gradient>1/6;
figure(7)
fill(Ex',Ey',real(Gradkrit)')
axis ('equal')

%% Med lerlager

clear all
close all
format compact
format short

b=[20.5];

%Punkter
Vertices = [
    0 0 %1
    10 0 %2
    20 0 %3
    20+b 0 %4
    30+b 0 %5
    0 3 %6
    10 3 %7
    20 3 %8
    20+b 3 %9
    30+b 3 %10
    0 6 %11
    10 6 %12
    20 6 %13
    20+b 6 %14
    30+b 6 %15
    0 10 %16
    10 10 %17
    20 10 %18
    20+b-0.025 10 %19
    20+b+0.025 10 %20
    30+b 10 %21
    ];

%Linjer
Segments = [
    1 2 %1
    2 3 %2
    3 4 %3
    4 5 %4
    6 7 %5
    7 8 %6
    8 9 %7
    9 10 %8
    11 12 %9
    12 13 %10
    13 14 %11
    14 15 %12
    16 17 %13
    17 18 %14
    18 19 %15
    20 21 %16
    1 6 %17
    6 11 %18
    11 16 %19
    2 7 %20
    7 12 %21
    12 17 %22
    3 8 %23
    8 13 %24
    13 18 %25
    4 9 %36
    9 14 %27
    14 19 %28
    14 20 %29
    5 10 %30
    10 15 %31
    15 21 %32
    ];

%Ytor
Surfaces=[
    1 20 5 17 %1
    2 23 6 20 %2
    3 26 7 23 %3
    4 30 8 26 %4
    5 21 9 18 %5
    6 24 10 21 %6
    7 27 11 24 %7
    8 31 12 27 %8
    9 22 13 19 %9
    10 25 14 22 %10
    11 28 15 25 %11
    12 32 16 29 %12
    ];

%Antal indelningar per linje
Seed =round([10 10 2*b 20 10 10 2*b 20 10 10 2*b 20 10 10 2*b 20 3 6 8 3 6 8 3 6 8 3 6 8 8 3 6 8]);

segp=[Seed];
nen=4;
dofsPerNode=1;
mp=[dofsPerNode, nen];

%Illustrerar elementindelning
geomdraw2(Vertices,Segments,Surfaces,segp,mp)

%Generera mesh
[Coord Edof Dof meshdb ]=strMeshgen(Vertices,Segments,Surfaces,segp,mp);

%Generera elementkoordinater
[Ex,Ey]=coordxtr(Edof,Coord,Dof,nen);

%Rita mesh med numrerade element
figure(2)
eldraw2(Ex,Ey,[1 4 1],Edof(:,1));

%Rita mesh med numrerade frihetsgrader
figure(3)
eldraw2(Ex,Ey,[1 4 0]);
text(Coord(:,1),Coord(:,2),num2str(Dof(:)))

%Skapa elementstyvhet, assemblera global styvhetsmatris, bestäm randvilkor
%och lös ekvationssystemet

%Egenskaper
ep=[1];
D=eye(2)*10^-4;
nel=size(Edof,1);
ndof=max(Edof(:));

%Nollor tas bort, elementstyvhetsmatris och assemblering
K=sparse(ndof,ndof);
f=sparse(ndof,1);
for i=1:nel
    Ke=flw2qe(Ex(i,:),Ey(i,:),ep,D);
    K=sparse_assem(Edof(i,:),K,Ke);
end

%Randvillkor
bc1=extrSeg([13 14]',meshdb,1);
bc2=extrSeg([16]',meshdb,1);

%Definera randvillkor
bc=[bc1, ones(length(bc1),1)*7; bc2, ones(length(bc2),1)*0];

%Lös ekvationssystemet
[a,r]=solveq(K,f,bc);
Ed=extract_ed(Edof,a);

%Beräkna element flux och gradient
for i=1:nel
    [Es(i,:),Et(i,:)]=flw2qs(Ex(i,:),Ey(i,:),ep,D,Ed(i,:));
end

%Plotta flödespilar
figure(4)
eldraw2(Ex,Ey,[1 4 0]);
[sfac]=elflux2(Ex(550,:),Ey(550,:),Es(550,:));
elflux2(Ex,Ey,Es,[1 2],sfac);

%Plotta tryck
figure(5)
eldraw2(Ex,Ey, [1 4 0]);
eliso2(Ex,Ey,Ed,10,[1 2 1]);

%Plotta färger
figure(6)
fill(Ex',Ey',Ed')
axis('equal')

%Flöde = 1.9825e-04
flodein=sum(r(bc1))
flodeut=sum(r(bc2))

%Gradient risk för erosion
[Gradient]=((Et(:,2).^2+Et(:,1).^2).^0.5);
Gradkrit=Gradient>1/6;
figure(7)
fill(Ex',Ey',real(Gradkrit)')
axis ('equal')