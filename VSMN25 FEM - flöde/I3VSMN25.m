% Assignment 3 - Transient Heat Flux - Isabella McAuley Skriver & Andrea Nilsson
%Fel i rad 195 och 479

%% Del a

close all
format compact

% Konvektionslager
h=0.005;

% Tidssteg, dt=1 timme
% Tidsintegrationskonstant, (0, 0.5, 1)
alpha=1;

% Definiera vertexkoordinater
Vertices=[0 0;                                                  % 1
    0.075*cos(30*pi/180) 0                                      % 2 
    0.210+0.075*cos(150*pi/180) 0;                              % 3 
    0.210 0;                                                    % 4
    0 0.035;                                                    % 5
    0.210 0.035;                                                % 6
    0.075*cos(330*pi/180), 0.110+0.075*sin(330*pi/180);         % 7 
    0.210+0.075*cos(210*pi/180), 0.110+0.075*sin(210*pi/180);   % 8  
    0.075*cos(30*pi/180), 0.110+0.075*sin(30*pi/180);           % 9 
    0.210+0.075*cos(150*pi/180), 0.110+0.075*sin(150*pi/180);   % 10
    0 0.185;                                                    % 11
    0.210 0.185;                                                % 12
    0 0.220;                                                    % 13
    0.075*cos(30*pi/180), 0.220;                                % 14
    0.210+0.075*cos(150*pi/180) 0.220;                          % 15
    0.210 0.220;                                                % 16
    0 0.220+h;                                                  % 17
    0.075*cos(30*pi/180) 0.220+h;                               % 18
    0.210+0.075*cos(150*pi/180) 0.220+h;                        % 19
    0.210 0.220+h;                                              % 20
    0.075*cos(300*pi/180), 0.110+0.075*sin(300*pi/180);         % 21 extra
    0.075*cos(0*pi/180), 0.110+0.075*sin(0*pi/180);             % 22 extra
    0.075*cos(60*pi/180), 0.110+0.075*sin(60*pi/180);           % 23 extra
    0.210+0.075*cos(240*pi/180), 0.110+0.075*sin(240*pi/180);   % 24 extra
    0.210+0.075*cos(180*pi/180), 0.110+0.075*sin(180*pi/180);   % 25 extra
    0.210+0.075*cos(120*pi/180), 0.110+0.075*sin(120*pi/180);]; % 26 extra

% Definiera linjesegment
Segments=[1 2;      % 1
    2 3;            % 2
    3 4;            % 3
    1 5;            % 4
    2 7;            % 5
    3 8;            % 6
    4 6;            % 7
    5 7;            % 8
    7 8;            % 9
    8 6;            % 10
    7 9;            % 11
    8 10;           % 12
    6 12;           % 13
    11 9;           % 14
    9 10;           % 15
    10 12;          % 16
    11 13;          % 17
    9 14;           % 18
    10 15;          % 19
    12 16;          % 20
    13 14;          % 21
    14 15;          % 22
    15 16;          % 23
    13 17;          % 24
    14 18;          % 25
    15 19;          % 26
    16 20;          % 27
    17 18;          % 28
    18 19;          % 29
    19 20];         % 30

% Definiera ytor utifrån segmentnummer
Surfaces=[1 5 8 4;  % 1
    2 6 9 5;        % 2
    3 7 10 6;       % 3
    9 12 15 11;     % 4
    10 13 16 12;    % 5
    14 18 21 17;    % 6
    15 19 22 18;    % 7
    16 20 23 19;    % 8
    21 25 28 24;    % 9
    22 26 29 25;    % 10
    23 27 30 26];   % 11

% Definiera antalet element på varje segment
Seed=round([7 10 7 6 6 6 6 7 10 7 7 7 7 7 10 7 6 6 6 6 7 10 7 1 1 1 1 7 10 7]);

% Definiera kurvorna
iso8=zeros(1,30); 
iso8(8)=21;     iso8(11)=22;    iso8(14)=23;
iso8(10)=24;    iso8(12)=25;    iso8(16)=26;

% Kombinera seed och kurvor
Segp(1:2:60)=Seed;  Segp(2:2:60)=iso8;

nen=3;  % Triangelelement
dofsPerNode=1; % 1 dof per nod
mp=[dofsPerNode, nen];

% Rita geometrin
geomdraw2(Vertices, Segments, Surfaces, Segp, mp)

% Skapa mesh
[Coord Edof Dof meshdb]=strMeshgen(Vertices, Segments, Surfaces, Segp, mp);

% Skapa elementkoordinater
[Ex, Ey]=coordxtr(Edof, Coord, Dof, nen);

% Rita element mesh med numrering av elementen
figure(2)
eldraw2(Ex, Ey, [1 4 0])

% Skapa element styvhetsmatris, assemblera till global styvhetsmatris,
% skapa randvillkor och lös ekvationssystemet
nel=size(Edof,1);
ndof=length(Dof);
%
K=sparse(ndof,ndof);
C=sparse(ndof,ndof);
f=sparse(ndof,1);

% Skapa styvhetsmatris för betongelement och assemblera till global
% styvhetsmatris
for j=[1 2 3 4 6 7 8] % Ytor med betong
    t=1;
    rho=2400;
    c=1000; 
    ep=[t, rho, c];
    k=1.7;
    Dc=k*[1 0; 0 1];
    elems=extrSurf(j,meshdb);
    for k=elems'
        [Ke, Ce]=flw2tt(Ex(k,:),Ey(k,:),ep,Dc);
        K=sparse_assem(Edof(k,:),K,Ke);
        C=sparse_assem(Edof(k,:),C,Ce);
    end
end

% Skapa styvhetsmatris för element med hålrum och assemblera till 
% global styvhetsmatris
for j=5     % Ytor med luft
    t=1;
    rho=1.29;
    c=1000;
    ep=[t, rho,c];
    k=0.026;
    Da=k*[1 0; 0 1];
    elems=extrSurf(j,meshdb);
    for k=elems'
        [Ke, Ce]=flw2tt(Ex(k,:),Ey(k,:),ep,Da);
        K=sparse_assem(Edof(k,:),K,Ke);
        C=sparse_assem(Edof(k,:),C,Ce);
    end
end

% Skapa styvhetsmatris för element med konvektion och assemblera till 
% global styvhetsmatris
for j=[9 10 11] % Ytor med konvektionslager
    t=1;
    rho=1.29;
    c=1000;
    ep=[t, rho, c];
    alphac=1/0.13;
    Dk=alphac*h*[0 0; 0 1];        % h är tjocklek på fiktivt konvektionslager
    elems=extrSurf(j,meshdb);
    for k=elems'
        [Ke, Ce]=flw2tt(Ex(k,:),Ey(k,:),ep,Dk);
        K=sparse_assem(Edof(k,:),K,Ke);
        C=sparse_assem(Edof(k,:),C,Ce);
    end
end

% Randvillkor, 
% innetemperatur 21 grader, temperatur i hålrum 28 grader
bc1=extrSeg([28, 29, 30]',meshdb,1);
bc2=extrSeg([8, 11, 14]',Segments,Segp,Vertices,Dof,Coord,1);

bc=[bc1, ones(length(bc1),1)*21; bc2, ones(length(bc2),1)*28];

% Definiera initiell temperatur vid t=0
a0=ones(ndof,1)*10;

% Definiera och utför tidsstegring under 16 timmar (sekunder)
% Definiera indata för tidsstegringskommandot
ip=[3600 3600*16 alpha];

% Tider då resultatet extraheras
snap=[3600:16*3600];

% Utför time-stepping
[a,da]=step1(K,C,f,a0,bc,ip);

% Extrahera temperaturer efter 1, 4, 8 och 16 timmar
Ed1=extract_ed(Edof,a(:,2));
Ed4=extract_ed(Edof,a(:,5));
Ed8=extract_ed(Edof,a(:,9));
Ed16=extract_ed(Edof,a(:,17));

% Beräkna heat flux 
for j=1:11
    if ((j==9) | (j==10) | (j==11))
        Dm=Dk;
    elseif ((j==5))
        Dm=Da;
    else
        Dm=Dc;
    end
    elems=extrSurf(j,meshdb);
    for i=elems'
        [Es1(i,:), Et1(i,:)]=flw2ts(Ex(i,:), Ey(i,:),Dm,Ed1(i,:));
        [Es4(i,:),Et4(i,:)]=flw2ts(Ex(i,:),Ey(i,:),Dm,Ed4(i,:));
        [Es8(i,:),Et8(i,:)]=flw2ts(Ex(i,:),Ey(i,:),Dm,Ed8(i,:));
        [Es16(i,:),Et16(i,:)]=flw2ts(Ex(i,:),Ey(i,:),Dm,Ed16(i,:));
    end
end

figure(4)
subplot(2,2,1)
fill(Ex',Ey',Ed1')
axis("equal")
title('Temperaturfördelning efter 1 timme')

subplot(2,2,2)
fill(Ex',Ey',Ed4')
axis("equal")
title('Temperaturfördelning efter 4 timmar')

subplot(2,2,3)
fill(Ex',Ey',Ed8')
axis("equal")
title('Temperaturfördelning efter 8 timmar')

subplot(2,2,4)
fill(Ex',Ey',Ed16')
axis("equal")
title('Temperaturfördelning efter 16 timmar')

% Temperaturfördelning längs betongens övre yta
[top_dof top_vertex]=extrSeg([21 22 23]',meshdb,1);
koord=Coord(top_vertex,1);
[koord_sort,koord_index]=sort(koord);
T1=a(top_vertex,2);
T1=T1(koord_index);
T4=a(top_vertex,5);
T4=T4(koord_index);
T8=a(top_vertex,9);
T8=T8(koord_index);
T16=a(top_vertex,17);
T16=T16(koord_index);
figure(5)
plot(koord_sort,T1,'b',koord_sort,T4,'r',koord_sort,T8,'g',koord_sort,T16,'k')
legend('Temperatur efter 1 timme','Temperatur efter 4 timmar','Temperatur efter 8 timmar','Temperatur efter 16 timmar')
grid on
axis([0 0.21 10 30])
xlabel('Längd längs betongens övre yta [m]')
ylabel('Temperatur [°C]')
title('Temperaturfördelning längs betongens övre yta')

% Beräknar transient värmeflöde in i rummet per kvadratmeter golvyta som en
% funktion av tiden
f=K*a+C*da;
[top_dof top_vertex]=extrSeg([28 29 30]',meshdb,1);
frum=f(top_dof,:);
heatflow=sum(frum);
heatflow=heatflow/0.21;

% Plottar total heat flow into the room per m^2
figure(6)
timmar=0:16;
plot(timmar,heatflow)
xlabel('Timmar')
ylabel('Värmeflöde per area')
title('Värmeflöde per area mellan betongen och rummet')

% Hur det totala flödet varierar med tiden. (steady-state)
flode=sum(f);
figure(7)
plot(timmar,flode)
xlabel('Timmar')
ylabel('Totalt flöde')
title('Värmeflödets variation över tid')


%% Del b

close all
format compact

% Konvektionslager
h=0.005;

% Tidssteg, dt=60s
% Tidsintegrationskonstant, (0, 0.5, 1)
alpha=0.5;

% Definiera vertexkoordinater
Vertices=[0 0;                                                  % 1
    0.075*cos(30*pi/180) 0                                      % 2 
    0.210+0.075*cos(150*pi/180) 0;                              % 3 
    0.210 0;                                                    % 4
    0 0.035;                                                    % 5
    0.210 0.035;                                                % 6
    0.075*cos(330*pi/180), 0.110+0.075*sin(330*pi/180);         % 7 
    0.210+0.075*cos(210*pi/180), 0.110+0.075*sin(210*pi/180);   % 8  
    0.075*cos(30*pi/180), 0.110+0.075*sin(30*pi/180);           % 9 
    0.210+0.075*cos(150*pi/180), 0.110+0.075*sin(150*pi/180);   % 10
    0 0.185;                                                    % 11
    0.210 0.185;                                                % 12
    0 0.220;                                                    % 13
    0.075*cos(30*pi/180), 0.220;                                % 14
    0.210+0.075*cos(150*pi/180) 0.220;                          % 15
    0.210 0.220;                                                % 16
    0 0.220+h;                                                  % 17
    0.075*cos(30*pi/180) 0.220+h;                               % 18
    0.210+0.075*cos(150*pi/180) 0.220+h;                        % 19
    0.210 0.220+h;                                              % 20
    0.075*cos(300*pi/180), 0.110+0.075*sin(300*pi/180);         % 21 extra
    0.075*cos(0*pi/180), 0.110+0.075*sin(0*pi/180);             % 22 extra
    0.075*cos(60*pi/180), 0.110+0.075*sin(60*pi/180);           % 23 extra
    0.210+0.075*cos(240*pi/180), 0.110+0.075*sin(240*pi/180);   % 24 extra
    0.210+0.075*cos(180*pi/180), 0.110+0.075*sin(180*pi/180);   % 25 extra
    0.210+0.075*cos(120*pi/180), 0.110+0.075*sin(120*pi/180);]; % 26 extra

% Definiera linjesegment
Segments=[1 2;      % 1
    2 3;            % 2
    3 4;            % 3
    1 5;            % 4
    2 7;            % 5
    3 8;            % 6
    4 6;            % 7
    5 7;            % 8
    7 8;            % 9
    8 6;            % 10
    7 9;            % 11
    8 10;           % 12
    6 12;           % 13
    11 9;           % 14
    9 10;           % 15
    10 12;          % 16
    11 13;          % 17
    9 14;           % 18
    10 15;          % 19
    12 16;          % 20
    13 14;          % 21
    14 15;          % 22
    15 16;          % 23
    13 17;          % 24
    14 18;          % 25
    15 19;          % 26
    16 20;          % 27
    17 18;          % 28
    18 19;          % 29
    19 20];         % 30

% Definiera ytor utifrån segmentnummer
Surfaces=[1 5 8 4;  % 1
    2 6 9 5;        % 2
    3 7 10 6;       % 3
    9 12 15 11;     % 4
    10 13 16 12;    % 5
    14 18 21 17;    % 6
    15 19 22 18;    % 7
    16 20 23 19;    % 8
    21 25 28 24;    % 9
    22 26 29 25;    % 10
    23 27 30 26];   % 11

% Definiera antalet element på varje segment
Seed=round([7 10 7 6 6 6 6 7 10 7 7 7 7 7 10 7 6 6 6 6 7 10 7 1 1 1 1 7 10 7]);

% Definiera kurvorna
iso8=zeros(1,30); 
iso8(8)=21;     iso8(11)=22;    iso8(14)=23;
iso8(10)=24;    iso8(12)=25;    iso8(16)=26;

% Kombinera seed och kurvor
Segp(1:2:60)=Seed;  Segp(2:2:60)=iso8;

nen=3;  % Triangelelement
dofsPerNode=1; % 1 dof per nod
mp=[dofsPerNode, nen];

% Rita geometrin
geomdraw2(Vertices, Segments, Surfaces, Segp, mp)

% Skapa mesh
[Coord Edof Dof meshdb]=strMeshgen(Vertices, Segments, Surfaces, Segp, mp);

% Skapa elementkoordinater
[Ex, Ey]=coordxtr(Edof, Coord, Dof, nen);

% Rita element mesh med numrering av elementen
figure(2)
eldraw2(Ex, Ey, [1 4 0])

% Skapa element styvhetsmatris, assemblera till global styvhetsmatris,
% skapa randvillkor och lös ekvationssystemet
nel=size(Edof,1);
ndof=length(Dof);
%
K=sparse(ndof,ndof);
C=sparse(ndof,ndof);
f=sparse(ndof,1);

% Skapa styvhetsmatris för betongelement och assemblera till global
% styvhetsmatris
for j=[1 2 3 4 6 7 8] % Ytor med betong
    t=1;
    rho=2400;
    c=1000; 
    ep=[t, rho, c];
    k=1.7;
    Dc=k*[1 0; 0 1];
    elems=extrSurf(j,meshdb);
    for k=elems'
        [Ke, Ce]=flw2tt(Ex(k,:),Ey(k,:),ep,Dc);
        K=sparse_assem(Edof(k,:),K,Ke);
        C=sparse_assem(Edof(k,:),C,Ce);
    end
end

% Skapa styvhetsmatris för element med hålrum och assemblera till 
% global styvhetsmatris
for j=5     % Ytor med luft
    t=1;
    rho=1.29;
    c=1000;
    ep=[t, rho,c];
    k=0.026;
    Da=k*[1 0; 0 1];
    elems=extrSurf(j,meshdb);
    for k=elems'
        [Ke, Ce]=flw2tt(Ex(k,:),Ey(k,:),ep,Da);
        K=sparse_assem(Edof(k,:),K,Ke);
        C=sparse_assem(Edof(k,:),C,Ce);
    end
end

% Skapa styvhetsmatris för element med konvektion och assemblera till 
% global styvhetsmatris
for j=[9 10 11] % Ytor med konvektionslager
    t=1;
    rho=1.29;
    c=1000;
    ep=[t, rho, c];
    alphac=1/0.13;
    Dk=alphac*h*[0 0; 0 1];        % h är tjocklek på fiktivt konvektionslager
    elems=extrSurf(j,meshdb);
    for k=elems'
        [Ke, Ce]=flw2tt(Ex(k,:),Ey(k,:),ep,Dk);
        K=sparse_assem(Edof(k,:),K,Ke);
        C=sparse_assem(Edof(k,:),C,Ce);
    end
end

% Randvillkor, 
% innetemperatur 21 grader, temperatur i hålrum 28 grader
bc1=extrSeg([28, 29, 30]',meshdb,1);
bc2=extrSeg([8, 11, 14]',Segments,Segp,Vertices,Dof,Coord,1);

bc=[bc1, ones(length(bc1),1)*21; bc2, ones(length(bc2),1)*28];

% Definiera initiell temperatur vid t=0
a0=ones(ndof,1)*10;

% Definiera och utför tidsstegring under 16 timmar (sekunder)
% Definiera indata för tidsstegringskommandot
ip=[60 3600*16 alpha];

% Tider då resultatet extraheras
snap=[3600:16*3600];

% Utför time-stepping
[a,da]=step1(K,C,f,a0,bc,ip);

% Extrahera temperaturer efter 1, 4, 8 och 16 timmar
Ed1=extract_ed(Edof,a(:,61));
Ed4=extract_ed(Edof,a(:,241));
Ed8=extract_ed(Edof,a(:,481));
Ed16=extract_ed(Edof,a(:,961));

% Beräkna heat flux 
for j=1:11
    if ((j==9) | (j==10) | (j==11))
        Dm=Dk;
    elseif ((j==5))
        Dm=Da;
    else
        Dm=Dc;
    end
    elems=extrSurf(j,meshdb);
    for i=elems'
        [Es1(i,:), Et1(i,:)]=flw2ts(Ex(i,:), Ey(i,:),Dm,Ed1(i,:));
        [Es4(i,:),Et4(i,:)]=flw2ts(Ex(i,:),Ey(i,:),Dm,Ed4(i,:));
        [Es8(i,:),Et8(i,:)]=flw2ts(Ex(i,:),Ey(i,:),Dm,Ed8(i,:));
        [Es16(i,:),Et16(i,:)]=flw2ts(Ex(i,:),Ey(i,:),Dm,Ed16(i,:));
    end
end

figure(4)
subplot(2,2,1)
fill(Ex',Ey',Ed1')
axis("equal")
title('Temperaturfördelning efter 1 timme')

subplot(2,2,2)
fill(Ex',Ey',Ed4')
axis("equal")
title('Temperaturfördelning efter 4 timmar')

subplot(2,2,3)
fill(Ex',Ey',Ed8')
axis("equal")
title('Temperaturfördelning efter 8 timmar')

subplot(2,2,4)
fill(Ex',Ey',Ed16')
axis("equal")
title('Temperaturfördelning efter 16 timmar')

% Temperaturfördelning längs betongens övre yta
[top_dof top_vertex]=extrSeg([21 22 23]',meshdb,1);
koord=Coord(top_vertex,1);
[koord_sort,koord_index]=sort(koord);
T1=a(top_vertex,61);
T1=T1(koord_index);
T4=a(top_vertex,241);
T4=T4(koord_index);
T8=a(top_vertex,481);
T8=T8(koord_index);
T16=a(top_vertex,961);
T16=T16(koord_index);
figure(5)
plot(koord_sort,T1,'b',koord_sort,T4,'r',koord_sort,T8,'g',koord_sort,T16,'k')
legend('Temperatur efter 1 timme','Temperatur efter 4 timmar','Temperatur efter 8 timmar','Temperatur efter 16 timmar')
grid on
axis([0 0.21 12 32])
xlabel('Längd längs betongens övre yta [m]')
ylabel('Temperatur [°C]')
title('Temperaturfördelning längs betongens övre yta')

% Beräknar transient värmeflöde in i rummet per kvadratmeter golvyta som en
% funktion av tiden
f=K*a+C*da;
[top_dof top_vertex]=extrSeg([28 29 30]',meshdb,1);
frum=f(top_dof,:);
heatflow=sum(frum);
heatflow=heatflow/0.21;
% Plottar total heat flow into the room per m^2
figure(6)
minuter=0:960
plot(minuter,heatflow)
xlabel('Minuter')
ylabel('Värmeflöde per area golv')
title('Värmeflöde per area mellan betong och rum')

% Hur det totala flödet varierar med tiden. (steady-state)
flode=sum(f);
figure(7)
plot(minuter,flode)
xlabel('Minuter')
ylabel('Totalt flöde')
title('Värmeflödets variation över tid')

%%
close all
format compact

% Konvektionslager
h=0.005;

% Tidssteg, dt=1s
% Tidsintegrationskonstant, (0, 0.5, 1)
alpha=0.5;

% Definiera vertexkoordinater
Vertices=[0 0;                                                  % 1
    0.075*cos(30*pi/180) 0                                      % 2 
    0.210+0.075*cos(150*pi/180) 0;                              % 3 
    0.210 0;                                                    % 4
    0 0.035;                                                    % 5
    0.210 0.035;                                                % 6
    0.075*cos(330*pi/180), 0.110+0.075*sin(330*pi/180);         % 7 
    0.210+0.075*cos(210*pi/180), 0.110+0.075*sin(210*pi/180);   % 8  
    0.075*cos(30*pi/180), 0.110+0.075*sin(30*pi/180);           % 9 
    0.210+0.075*cos(150*pi/180), 0.110+0.075*sin(150*pi/180);   % 10
    0 0.185;                                                    % 11
    0.210 0.185;                                                % 12
    0 0.220;                                                    % 13
    0.075*cos(30*pi/180), 0.220;                                % 14
    0.210+0.075*cos(150*pi/180) 0.220;                          % 15
    0.210 0.220;                                                % 16
    0 0.220+h;                                                  % 17
    0.075*cos(30*pi/180) 0.220+h;                               % 18
    0.210+0.075*cos(150*pi/180) 0.220+h;                        % 19
    0.210 0.220+h;                                              % 20
    0.075*cos(300*pi/180), 0.110+0.075*sin(300*pi/180);         % 21 extra
    0.075*cos(0*pi/180), 0.110+0.075*sin(0*pi/180);             % 22 extra
    0.075*cos(60*pi/180), 0.110+0.075*sin(60*pi/180);           % 23 extra
    0.210+0.075*cos(240*pi/180), 0.110+0.075*sin(240*pi/180);   % 24 extra
    0.210+0.075*cos(180*pi/180), 0.110+0.075*sin(180*pi/180);   % 25 extra
    0.210+0.075*cos(120*pi/180), 0.110+0.075*sin(120*pi/180);]; % 26 extra

% Definiera linjesegment
Segments=[1 2;      % 1
    2 3;            % 2
    3 4;            % 3
    1 5;            % 4
    2 7;            % 5
    3 8;            % 6
    4 6;            % 7
    5 7;            % 8
    7 8;            % 9
    8 6;            % 10
    7 9;            % 11
    8 10;           % 12
    6 12;           % 13
    11 9;           % 14
    9 10;           % 15
    10 12;          % 16
    11 13;          % 17
    9 14;           % 18
    10 15;          % 19
    12 16;          % 20
    13 14;          % 21
    14 15;          % 22
    15 16;          % 23
    13 17;          % 24
    14 18;          % 25
    15 19;          % 26
    16 20;          % 27
    17 18;          % 28
    18 19;          % 29
    19 20];         % 30

% Definiera ytor utifrån segmentnummer
Surfaces=[1 5 8 4;  % 1
    2 6 9 5;        % 2
    3 7 10 6;       % 3
    9 12 15 11;     % 4
    10 13 16 12;    % 5
    14 18 21 17;    % 6
    15 19 22 18;    % 7
    16 20 23 19;    % 8
    21 25 28 24;    % 9
    22 26 29 25;    % 10
    23 27 30 26];   % 11

% Definiera antalet element på varje segment
Seed=round([7 10 7 6 6 6 6 7 10 7 7 7 7 7 10 7 6 6 6 6 7 10 7 1 1 1 1 7 10 7]);

% Definiera kurvorna
iso8=zeros(1,30); 
iso8(8)=21;     iso8(11)=22;    iso8(14)=23;
iso8(10)=24;    iso8(12)=25;    iso8(16)=26;

% Kombinera seed och kurvor
Segp(1:2:60)=Seed;  Segp(2:2:60)=iso8;

nen=3;  % Triangelelement
dofsPerNode=1; % 1 dof per nod
mp=[dofsPerNode, nen];

% Rita geometrin
geomdraw2(Vertices, Segments, Surfaces, Segp, mp)

% Skapa mesh
[Coord Edof Dof meshdb]=strMeshgen(Vertices, Segments, Surfaces, Segp, mp);

% Skapa elementkoordinater
[Ex, Ey]=coordxtr(Edof, Coord, Dof, nen);

% Rita element mesh med numrering av elementen
figure(2)
eldraw2(Ex, Ey, [1 4 0])

% Skapa element styvhetsmatris, assemblera till global styvhetsmatris,
% skapa randvillkor och lös ekvationssystemet
nel=size(Edof,1);
ndof=length(Dof);

%
K=sparse(ndof,ndof);
C=sparse(ndof,ndof);
f=sparse(ndof,1);

% Skapa styvhetsmatris för betongelement och assemblera till global
% styvhetsmatris
for j=[1 2 3 4 6 7 8] % Ytor med betong
    t=1;
    rho=2400;
    c=1000; 
    ep=[t, rho, c];
    k=1.7;
    Dc=k*[1 0; 0 1];
    elems=extrSurf(j,meshdb);
    for k=elems'
        [Ke, Ce]=flw2tt(Ex(k,:),Ey(k,:),ep,Dc);
        K=sparse_assem(Edof(k,:),K,Ke);
        C=sparse_assem(Edof(k,:),C,Ce);
    end
end

% Skapa styvhetsmatris för element med hålrum och assemblera till 
% global styvhetsmatris
for j=5     % Ytor med luft
    t=1;
    rho=1.29;
    c=1000;
    ep=[t, rho,c];
    k=0.026;
    Da=k*[1 0; 0 1];
    elems=extrSurf(j,meshdb);
    for k=elems'
        [Ke, Ce]=flw2tt(Ex(k,:),Ey(k,:),ep,Da);
        K=sparse_assem(Edof(k,:),K,Ke);
        C=sparse_assem(Edof(k,:),C,Ce);
    end
end

% Skapa styvhetsmatris för element med konvektion och assemblera till 
% global styvhetsmatris
for j=[9 10 11] % Ytor med konvektionslager
    t=1;
    rho=1.29;
    c=1000;
    ep=[t, rho, c];
    alphac=1/0.13;
    Dk=alphac*h*[0 0; 0 1];        % h är tjocklek på fiktivt konvektionslager
    elems=extrSurf(j,meshdb);
    for k=elems'
        [Ke, Ce]=flw2tt(Ex(k,:),Ey(k,:),ep,Dk);
        K=sparse_assem(Edof(k,:),K,Ke);
        C=sparse_assem(Edof(k,:),C,Ce);
    end
end

% Randvillkor, 
% innetemperatur 21 grader, temperatur i hålrum 28 grader
bc1=extrSeg([28, 29, 30]',meshdb,1);
bc2=extrSeg([8, 11, 14]',Segments,Segp,Vertices,Dof,Coord,1);

bc=[bc1, ones(length(bc1),1)*21; bc2, ones(length(bc2),1)*28];

% Definiera initiell temperatur vid t=0
a0=ones(ndof,1)*10;

% Definiera och utför tidsstegring under 16 timmar (sekunder)
% Definiera indata för tidsstegringskommandot
ip=[1 3600*16 alpha];

% Tider då resultatet extraheras
snap=[3600:16*3600];

% Utför time-stepping
[a,da]=step1(K,C,f,a0,bc,ip);

% Extrahera temperaturer efter 1, 4, 8 och 16 timmar
Ed1=extract_ed(Edof,a(:,3601));
Ed4=extract_ed(Edof,a(:,14401));
Ed8=extract_ed(Edof,a(:,28801));
Ed16=extract_ed(Edof,a(:,57601));

% Beräkna heat flux 
for j=1:11
    if ((j==9) | (j==10) | (j==11))
        Dm=Dk;
    elseif ((j==5))
        Dm=Da;
    else
        Dm=Dc;
    end
    elems=extrSurf(j,meshdb);
    for i=elems'
        [Es1(i,:), Et1(i,:)]=flw2ts(Ex(i,:), Ey(i,:),Dm,Ed1(i,:));
        [Es4(i,:),Et4(i,:)]=flw2ts(Ex(i,:),Ey(i,:),Dm,Ed4(i,:));
        [Es8(i,:),Et8(i,:)]=flw2ts(Ex(i,:),Ey(i,:),Dm,Ed8(i,:));
        [Es16(i,:),Et16(i,:)]=flw2ts(Ex(i,:),Ey(i,:),Dm,Ed16(i,:));
    end
end

figure(4)
subplot(2,2,1)
fill(Ex',Ey',Ed1')
axis("equal")
title('Temperaturfördelning efter 1 timme')

subplot(2,2,2)
fill(Ex',Ey',Ed4')
axis("equal")
title('Temperaturfördelning efter 4 timmar')

subplot(2,2,3)
fill(Ex',Ey',Ed8')
axis("equal")
title('Temperaturfördelning efter 8 timmar')

subplot(2,2,4)
fill(Ex',Ey',Ed16')
axis("equal")
title('Temperaturfördelning efter 16 timmar')

% Temperaturfördelning längs betongens övre yta
[top_dof top_vertex]=extrSeg([21 22 23]',meshdb,1);
koord=Coord(top_vertex,1);
[koord_sort,koord_index]=sort(koord);
T1=a(top_vertex,3601);
T1=T1(koord_index);
T4=a(top_vertex,14401);
T4=T4(koord_index);
T8=a(top_vertex,28801);
T8=T8(koord_index);
T16=a(top_vertex,57601);
T16=T16(koord_index);
figure(5)
plot(koord_sort,T1,'b',koord_sort,T4,'r',koord_sort,T8,'g',koord_sort,T16,'k')
legend('Temperatur efter 1 timme','Temperatur efter 4 timmar','Temperatur efter 8 timmar','Temperatur efter 16 timmar')
grid on
axis([0 0.21 10 30])
xlabel('Längd längs betongens övre yta [m]')
ylabel('Temperatur [°C]')
title('Temperaturfördelning längs betongens övre yta')

% Beräknar transient värmeflöde in i rummet per kvadratmeter golvyta som en
% funktion av tiden
f=K*a+C*da;
[top_dof top_vertex]=extrSeg([28 29 30]',meshdb,1);
frum=f(top_dof,:);
heatflow=sum(frum);
heatflow=heatflow/0.21;

% Plottar total heat flow into the room per m^2
figure(6)
sekunder=0:57600;
plot(sekunder,heatflow)
xlabel('Sekunder')
ylabel('Värmeflöde per area')
title('Värmeflöde per area mellan betongen och rummet')

% Hur det totala flödet varierar med tiden. (steady-state)
flode=sum(f);
figure(7)
plot(sekunder,flode)
axis([0 57600 0 100])
xlabel('Sekunder')
ylabel('Totalt flöde')
title('Värmeflödets variation över tid')

