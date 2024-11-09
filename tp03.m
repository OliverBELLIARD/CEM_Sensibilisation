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
V(25:28, 8:34) = v1; % Conducteur 1
V(5:22, 20:21) = v2; % Conducteur 2
figure(1)
subplot(1, 2, 1);
pcolor(V)
title("Potentiels avant itérations")

% Changer la palette de couleurs
colormap(jet)  % Palette allant du bleu au rouge
colorbar;  % Ajouter une barre de couleur

%% Calcul de convergence
Iter=0;     % Nombre d'itérations
seuil=1e-6; % Seuil de différence
cond=1;     % Condition de convergence
ii=2:Nx-1;
jj=2:Ny-1;

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
    V(25:28, 8:34) = v1; % Conducteur 1
    V(5:22, 20:21) = v2; % Conducteur 2

    % Equation de calcul
    V(ii,jj)=0.25*( V(ii+1,jj) + V(ii-1,jj) + V(ii,jj+1) + V(ii,jj-1) );

    % Calcul de condition de convergence
    cond=norm(abs(Vold(:)-V(:)));
    Iter=Iter+1;
end

% Temps écoulé
temps = toc;  % Fin du chronométrage

%% Figure
subplot(1, 2, 2);
pcolor(V)
title("Potentiels après "+Iter+" itérations", ...
    "Seuil : "+seuil+", Temps ecoulé : "+temps+"s")

% Changer la palette de couleurs
colormap(jet)  % Palette allant du bleu au rouge
colorbar;  % Ajouter une barre de couleur
