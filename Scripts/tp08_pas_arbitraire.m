% Résolution de l'équation de Laplace avec pas non unitaires
%
clear
close all
clc

%% Constantes
eps0 = 8.854e-12;  % Permittivité du vide en F/m

%% Dimensions / maillage
dx = 0.01;  % Pas en x (en cm)
dy = 0.01;  % Pas en y (en cm)

% Calcul du nombre de points en fonction de la taille du domaine et du pas
Lx = 40; % Largeur du domaine en cm
Ly = 40; % Hauteur du domaine en cm
Nx = round(Lx / dx); % Nombre de points en x
Ny = round(Ly / dy); % Nombre de points en y

%% Potentiels / sources
v0 = 0;   % Condition aux limites (en V)
v1 = 100; % Conducteur 1
v2 = -100; % Conducteur 2

% Initialisation de la matrice de calcul
V = zeros(Nx, Ny);  % Mettre toute la matrice à zéro

%% Sources
% Dimensions du potentiel 1
DimL1 = 28;
DimH1 = 4;
PotL1 = round(DimL1 / dy);
PotH1 = round(DimH1 / dy);
Y_abs_offset1 = 7;
Y_offset1 = round(Y_abs_offset1 / dy);

% Dimensions du potentiel 2
DimL2 = 2;
DimH2 = 18;
PotL2 = round(DimL2 / dy);
PotH2 = round(DimH2 / dy);
Y_abs_offset2 = -6;
Y_offset2 = round(Y_abs_offset2 / dy);

V(round((Y_offset1 + Ny / 2) - PotH1 / 2):round((Y_offset1 + Ny / 2) + PotH1 / 2 - 1), ...
    round((Nx / 2) - PotL1 / 2):round((Nx / 2) + PotL1 / 2)) = v1; % Conducteur 1 centré
V(round((Y_offset2 + Ny / 2) - PotH2 / 2):round((Y_offset2 + Ny / 2) + PotH2 / 2 - 1), ...
    round((Nx / 2) - PotL2 / 2):round((Nx / 2) + PotL2 / 2)) = v2; % Conducteur 2 centré

%% Calcul de convergence
Iter = 0;   % Nombre d'itérations
seuil = 1; % Seuil de différence
cond = 10;   % Condition de convergence
ii=1+1:Nx-1;
jj=1+1:Ny-1;

% Mesure du temps de calcul
tic;  % Début du chronométrage

while cond > seuil
    % Mémoire de la matrice précédente
    Vold = V;

    % Conditions aux limites
    V(1,:) = v0;    % Première colonne à zéro
    V(Nx,:) = v0;   % Dernière colonne à zéro
    V(:,1) = v0;    % Première ligne à zéro
    V(:,Ny) = v0;   % Dernière ligne à zéro

    % Sources
    V(round((Y_offset1 + Ny / 2) - PotH1 / 2):round((Y_offset1 + Ny / 2) + PotH1 / 2 - 1), ...
        round((Nx / 2) - PotL1 / 2):round((Nx / 2) + PotL1 / 2)) = v1; % Conducteur 1 centré
    V(round((Y_offset2 + Ny / 2) - PotH2 / 2):round((Y_offset2 + Ny / 2) + PotH2 / 2 - 1), ...
        round((Nx / 2) - PotL2 / 2):round((Nx / 2) + PotL2 / 2)) = v2; % Conducteur 2 centré

    % Equation de calcul avec des pas non unitaires
    V(ii, jj) = 0.25 * (V(ii+1, jj) + V(ii-1, jj) + V(ii, jj+1) + V(ii, jj-1));

    % Calcul de condition de convergence
    cond = norm(abs(Vold(:) - V(:)));
    Iter = Iter + 1;
end

%% Calcul du champ électrique
[Ex, Ey] = gradient(V, dx, dy);  % Calcul avec des pas dx et dy
Ex = -Ex; 
Ey = -Ey;

% Temps écoulé
temps = toc;  % Fin du chronométrage

%% Figure
figure(1)
% Figure du potentiel électrique
subplot(1, 2, 1);
h = pcolor(V);
set(h, 'EdgeColor', 'none'); % Hide grid

title("Potentiel V")
colormap(jet)  % Palette allant du bleu au rouge
colorbar;  % Ajouter une barre de couleur

% Visualisation des lignes de champ
subplot(1, 2, 2);
contour(V, 20);
% Nous pourrions afficher les vecteurs de champ mais elles nuisent à la
% visualisation du champ à cause de leur grand nombre.
%hold on
%quiver(Ex, Ey)

title("Champ V après "+Iter+" itérations", ...
    "Seuil : "+seuil+", Temps ecoulé : "+temps+"s")
colormap(jet)  % Palette allant du bleu au rouge
colorbar;  % Ajouter une barre de couleur

%% Calcul de la capacité
ContourY = round((Y_offset1 + Ny / 2) - PotH1 / 2 - 1):round((Y_offset1 + Ny / 2) + PotH1 / 2 + 1);
PosX_right = round(1 + (Nx / 2) + PotL1 / 2);
PosX_left = round(-1 + (Nx / 2) - PotL1 / 2);

% Calcul de l'intégrale de Gauss pour estimer Q1
Q1 = 0;

% Parcours du contour gauche
for x = ContourY
    Q1 = Q1 - eps0 * Ey(PosX_left, x) * dx;
end

% Parcours du contour droit
for x = ContourY
    Q1 = Q1 + eps0 * Ey(PosX_right, x) * dx;
end

% Parcours des côtés verticaux du contour
for y = PosX_left+1:PosX_right-1
    Q1 = Q1 - eps0 * Ex(y, ContourY(1)) * dy; % Côté gauche
    Q1 = Q1 + eps0 * Ex(y, ContourY(end)) * dy; % Côté droit
end

% Calcul de la capacité entre les deux conducteurs
C21 = Q1 / (v1 - v2);

% Affichage du résultat
disp(['La capacité Cij entre les deux conducteurs est : ', num2str(C21), ' F']);