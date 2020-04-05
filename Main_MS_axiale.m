%% Projet ST7 - Machine Axiale - ENE
%
%  objectif     :     Tracer la geometrie d'une machine synchrone axiale
%                     (type NS) à deux rotors et 1 stator
%                     avec du bobinage concentré
%
%  Auteur       :     Théo NDEREYIMANA, Lucas Roberto DAGORT, Matheus
%                     Elyasha LOPES, Valentin GOLDITÉ, Yoann ROUSSEL
%
%  Date de creation : 01-Avril-2020
%
%  Version      : v11

%--------------------------------------------------------------------------
%%                      Programmme Principal
%--------------------------------------------------------------------------

%% Adding Path to allow using femm through Matlab directly

clear all;
close all;
clc;
addpath('C:\femm42\mfiles') % Usual path

%% Discrétisation du modèle
pas = 20;                                                                  % Nombre de pas (1 pas du motor -> 2pi\(pas*P)) - P -> Paires de Pôles.

%% Declaration des paramètres

% Dimensions moteur

Lmax = 100;     %[mm]                                                      % Longueur de la machine.
Dext= 50;       %[mm]                                                      % Diametre exterieur.
Rext=Dext/2;    %[mm]                                                      % Rayon externe.
volume_encoche = 0.6;                                                      % 60% du volume du stator est destiné aux bobines.
volume_dents   = 0.4;                                                      % 40% du volume du stator est destiné aux dents.
Coeff_remplissage = 0.8;                                                   % Coeficient de remplissage des enconches.

% Alimentation et caracteristiques

Veff = 10;       % [Veff]                                                  % Tension d'alimentation (entrée de la machine pique-pique).
Ieff = 12;       % [Aeff]                                                  % Courant d'alimentation.
Imax_pique = 20; % [A]                                                     % Courant pique maximale.
%J = 8;          % [A/m^2]                                                 % Densité de courant.
N=10000;         % [rpm]                                                   % Vitesse de rotation désisé.

% Dimensions aimants

l_aimant = 2;          %[mm]                                               % longueur aimants.
d_aimant = 20;         %[mm]                                               % diametre aimants.
R_aimant = d_aimant/2; %[mm]                                               % Rayon aimant.

% Paramétres Matérieux

mu0=4*pi*1e-7;                                                             % parametres des materiaux.
Br0 = 1.16;                                                                % Induction remanente (20°C).
Hcb0 = 880*1e3;                                                            % Champ coercitive (20°C).
mua = Br0/(mu0*Hcb0);                                                      % Permeabilité relative de l'AP.
rho_cu = 1.72*1e-8;                                                        % Résisitivité du cuivre (ohm.m).
mvfer = 7800;                                                              % masse volumique du fer (kg/m3).
mvaim = 7500;                                                              % masse volumique du Aimants (kg/m3).
mvcu = 8960;                                                               % masse volumique du cuivre (kg/m3).
mvplastique = 1000;                                                        % masse volumique du plastique (kg/m3).
q = 1;                                                                     % Pertes fers spécifiques (W/kg).

%% Declaration des Variables

% Dimensions moteur

Dint = 10;      %[mm]                                                      % Diamètre interne.
Rint = Dint/2;  %[mm]                                                      % Rayon interne.
h_R = 15;       %[mm]                                                      % hauteur du rotor.
h_S = 19;       %[mm]                                                      % hauteur du stator.
h_e = 1;        %[mm]                                                      % hateur de l'entrefer.
h_b = 0.9*h_S;                                                             % Hauteur des bobines.
P = 2;                                                                     % Nombre de paires de poles.
nb_spires = 25;                                                            % Nombre de spires.

%% Linearization - 2D <- 3D  - pré-conditionnement

Req = (Rext+Rint)/2;                                                       % Rayon équivalent.
lm = 2*pi*Req;                                                             % Longuer équivalente.
prof_eq = (Rext^2-Rint^2)/(2*Req);                                         % Profondeur équivalente.

%% Entrées

Variables = [Rint, Req, lm, prof_eq, h_R, h_S, h_e, h_b, P, nb_spires];
Parameters = [Lmax, Rext, volume_encoche, volume_dents, Coeff_remplissage, Veff, Ieff,...
    Imax_pique, N, Hcb0, mua, rho_cu, mvfer, mvaim, mvcu, mvplastique, q, l_aimant, R_aimant];

%% Optimisation

lb = [4  5 15];                           % Lower bound
ub = [Rext-R_aimant 0.35*Lmax 0.28*Lmax]; % Upper bound
A = [];                                   % No linear inequality constraints
b = [];                                   % No linear inequality constraints
Aeq = [];                                 % No linear equality constraints
beq = [];                                 % No linear equality constraints

options = optimoptions(@gamultiobj,'PlotFcn',@gaplotpareto,'PopulationSize',5,'MaxGenerations',2);
[x,Fval,exitFlag,Output] = gamultiobj(@(x) myfunc(Parameters,pas,x),3,A,b,Aeq,beq,lb,ub,@(x) mycont(Parameters,pas,x),options);
myfunc(x)

function [Masse,Pertes] = myfunc(Parameters,pas,x)

%Paramètres
Lmax              = Parameters(1);
Rext              = Parameters(2);
volume_encoche    = Parameters(3);
volume_dents      = Parameters(4);
Coeff_remplissage = Parameters(5);
Ieff              = Parameters(7);
N                 = Parameters(9);
Hcb0              = Parameters(10);
mua               = Parameters(11);
rho_cu            = Parameters(12);
mvfer             = Parameters(13);
mvaim             = Parameters(14);                                                               
mvcu              = Parameters(15);                                                                 
mvplastique       = Parameters(16); 
q                 = Parameters(17);
l_aimant          = Parameters(18);
R_aimant          = Parameters(19);

% Variables
Rint      = x(1);
Req       = (Rext+Rint)/2;
lm        = 2*pi*Req;
prof_eq   = (Rext^2-Rint^2)/(2*Req);
h_R       = x(2);
h_S       = x(3);
h_e       = 1;
h_b       = 0.9*h_S;
P         = 2;
nb_spires = 25;

% Contraintes 
C1 = Rint + R_aimant - Rext;                                               % Rayon externe plus grand que la somme du rayin des aimants e du rayon interne.
C2 = l_aimant - h_R;                                                       % Hauteur du rotor plus grand que la hauteur des aimants.
C3 = 2*x(2) + h_S + 2*h_e - Lmax;                                          % La somme des hauteurs de chaque élement doit être plus petit que la longuer maximale.
C4 = 2*P*(2*R_aimant + 2) - 2*pi*Req;                                      % La distance entre les aimants doit être supérieur ou égal à 2 mm.

if C1<=0 && C2<=0 && C3<=0 && C4<=0
    Variables = [Rint, Req, lm, prof_eq, h_R, h_S, h_e, h_b, P, nb_spires];
    [~, ~,~, ~, Pertes, Masse] = MS_Axiale_Modele(Parameters, Variables, pas);
else
    disp('Machine infeasable')
end

end

function [c,ceq] = mycont(Parameters,pas,x)

% Paramétres
Lmax              = Parameters(1);
Rext              = Parameters(2);
volume_encoche    = Parameters(3);
volume_dents      = Parameters(4);
Coeff_remplissage = Parameters(5);
Ieff              = Parameters(7);
N                 = Parameters(9);
Hcb0              = Parameters(10);
mua               = Parameters(11);
rho_cu            = Parameters(12);
mvfer             = Parameters(13);
mvaim             = Parameters(14);                                                               
mvcu              = Parameters(15);                                                                 
mvplastique       = Parameters(16); 
q                 = Parameters(17);
l_aimant          = Parameters(18);
R_aimant          = Parameters(19);

% Variables
Rint      = x(1);
Req       = (Rext+Rint)/2;
lm        = 2*pi*Req;
prof_eq   = (Rext^2-Rint^2)/(2*Req);
h_R       = x(2);
h_S       = x(3);
h_e       = 1;
h_b       = 0.9*h_S;
P         = 2;
nb_spires = 25;

% Contraintes 
Cmin = 20;
FEM_max = 10;
J_max = 8;

C1 = Rint + R_aimant - Rext;                                               % Rayon externe plus grand que la somme du rayin des aimants e du rayon interne.
C2 = l_aimant - h_R;                                                       % Hauteur du rotor plus grand que la hauteur des aimants.
C3 = 2*h_R + h_S + 2*h_e - Lmax;                                           % La somme des hauteurs de chaque élement doit être plus petit que la longuer maximale.
C4 = 2*P*(2*R_aimant + 2) - 2*pi*Req;                                      % La distance entre les aimants doit être supérieur ou égal à 2 mm.

if C1<=0 && C2<=0 && C3<=0 && C4<=0
    
    Variables = [Rint, Req, lm, prof_eq, h_R, h_S, h_e, h_b, P, nb_spires];
    [~, ~, ~, Couple, ~, ~, FEM, J] = MS_Axiale_Modele(Parameters, Variables, pas);
    
    c(1,1) = Cmin - Couple;
    c(1,2) = FEM - FEM_max;
    c(1,3) = J - J_max;
    
else
    disp('Machine infeasable')
end
ceq = [];
end




