% Resolution equation de Laplace
%
clear
close all
clc

%% Dimensions / maillage
dx=1; % cm
dy=1; % cm
Nx = 40;
Ny = 40;

%% Potentiels / sources
v0 = 0; % condition aux limites (en V)
v1 = 100; % conducteur 1
v2 = -100; % conducteur 2

% Initialisation la matrice de calcul
V = zeros(Nx,Ny); % mettre toute la matrice a zero

%% Sources
% Dimensions du potentiel 1
Pot1L=28;
Pot1H=4;

% Dimensions du potentiel 2
Pot2L=2;
Pot2H=18;

V((7+Ny/2)-Pot1H/2:(7+Ny/2)+Pot1H/2-1, ...
    (Nx/2)-Pot1L/2:(Nx/2)+Pot1L/2) = v1; % Conducteur 1 centré
V((-6+Ny/2)-Pot2H/2:(-6+Ny/2)+Pot2H/2-1, ...
    (Nx/2)-Pot2L/2:(Nx/2)+Pot2L/2) = v2; % Conducteur 2 centré

%% Calcul de convergence
Iter=0;     % Nombre d'itérations
seuil=1e-2; % Seuil de différence
cond=1;     % Condition de convergence
ii=dx+1:Nx-dx;
jj=dy+1:Ny-dy;

% Mesure du temps de calcul
tic;  % Début du chronométrage

while cond>seuil
    % Mémoire de la matrice précédente
    Vold=V;

    % Conditions aux limites
    V(1,:) = v0;    % Première colone à zéro
    V(Nx,:) = v0;   % Dernière colone à zéro
    V(:,1) = v0;    % Première ligne à zéro
    V(:,Ny) = v0;   % Dernière ligne à zéro

    % Sources
    V((7+Ny/2)-Pot1H/2:(7+Ny/2)+Pot1H/2-1, ...
        (Nx/2)-Pot1L/2:(Nx/2)+Pot1L/2) = v1; % Conducteur 1 centré
    V((-6+Ny/2)-Pot2H/2:(-6+Ny/2)+Pot2H/2-1, ...
        (Nx/2)-Pot2L/2:(Nx/2)+Pot2L/2) = v2; % Conducteur 2 centré

    % Equation de calcul
    V(ii,jj)=0.25*( V(ii+1,jj) + V(ii-1,jj) + V(ii,jj+1) + V(ii,jj-1) );

    % Calcul de condition de convergence
    cond=norm(abs(Vold(:)-V(:)));
    Iter=Iter+1;
end

%% Calcul du champ
% Le champ électrique peut s'exprimer par la différence de deux potentiels
% en un point, c'est ce dont on se sert pour calculer nos champs (car dx et
% dy valent 1).
[Ex,Ey]=gradient(V);
Ex=-Ex; Ey=-Ey;

% Temps écoulé
temps = toc;  % Fin du chronométrage

%% Figure
figure(1)
subplot(1, 2, 1);
pcolor(V)
title("Potentiel V")

% Changer la palette de couleurs
colormap(jet)  % Palette allant du bleu au rouge
colorbar;  % Ajouter une barre de couleur

subplot(1, 2, 2);
contour(V)
hold on
quiver(Ex, Ey)
title("Champ V après "+Iter+" itérations", ...
    "Seuil : "+seuil+", Temps ecoulé : "+temps+"s")

% Changer la palette de couleurs
colormap(jet)  % Palette allant du bleu au rouge
colorbar;  % Ajouter une barre de couleur
