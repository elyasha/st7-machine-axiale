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
%  Version      : v10
%

%% Adding Path to allow using femm through Matlab directly

clear all; clc;

%addpath('C:\Program Files\femm42\mfiles') % Not the usual path
addpath('C:\femm42\mfiles') % Usual path


%% Open FEMM, Open New Document, Define the magnetostatic problem

openfemm();
newdocument(0);
mi_probdef(0,'millimeters','planar', 1E-8, 60, 30);

%% Declaration des parametres d'entrees

prompt = 'Selectionnez le type d´essaie: \n\n Essai à vide => 1  \n Essai avec charge => 2  \n \n Essai sélectionné: ';
Essai = input(prompt);

% Dimensions moteur
L = 100;        %[mm]                                                      % Longueur de la machine
Dint = 10;      %[mm]                                                      % Diamètre interne
Rint = Dint/2;  %[mm]                                                      % Rayon interne
Dext= 50;       %[mm]                                                      % Diametre exterieur
Rext=Dext/2;    %[mm]                                                      % Rayon externe
h_R = 15;       %[mm]                                                      % hauteur du rotor
h_S = 19;       %[mm]                                                      % hauteur du stator
h_e = 1;        %[mm]                                                      % hateur de l'entrefer

%--------------------------------------------------------------------------
P = 2;                                                                     % Nombre de paires de poles
nb_spires=25;                                                              % Nombre de spires
op=2*pi/(2*P);                                                             % Ouverture polaire

%--------------------------------------------------------------------------
% Dimensions aimants

l_aimant = 2;          %[mm]                                               % longueur aimants
d_aimant = 20;         %[mm]                                               % diametre aimants
R_aimant = d_aimant/2; % [mm]                                              % Rayon aimant

% Dimensions Bobines et dents

volume_encoche = 0.6;                                                      % 60% du volume entre les pôles est destiné aux bobines
volume_dents   = 0.4;                                                      % 40% du volume entre les pôles est destiné aux dents
h_b = 0.9*h_S  ;                                                           % Hauteur bobines
Coeff_remplissage = 0.8;                                                   % Coeficient de remplissage des enconches

%--------------------------------------------------------------------------
Surface_Rotor = (Rext^2-Rint^2)*pi;
Surface_aimant = 2*P*(pi*R_aimant^2);
Surf_rot_aimant = Surface_aimant/Surface_Rotor;

angle_aimant=((2*pi*Surf_rot_aimant)/(2*P))*(180/pi);                      % Surf_rot_aimant*100 % de la surface du rotor couverte par des aimants
angle_interaimant=((2*pi)*180/pi-angle_aimant*2*P)/(2*P);

%--------------------------------------------------------------------------
% Alimentation et caracteristiques

J = 8;           % [A/m^2]                                                 % Densité de courant
Veff = 10;       % [Veff]                                                  % Tension d'alimentation (entrée de la machine pique-pique)
Ieff = 10;       % [Aeff]
Imax_pique = 20; % [A]
psi = 0;         % [º]                                                     % Déphasage courant et tension d'entré
N=10000;         % [rpm]                                                   % Vitesse de rotation désisé
f=P*N/60;        % [Hz]                                                    % Frequence champ tournante
Tperiod=1/f;     % [s]                                                     % Periode champ tournante
w=2*pi*f;        % [rad/s]                                                 % Vitesse champ tournante
pas = 60;
%--------------------------------------------------------------------------
% Paramétres Matérieux

mu0=4*pi*1e-7;                                                             % parametres des materiaux
Br0 = 1.16;                                                                % Induction remanente (20°C)
Hcb0 = 880*1e3;                                                            % Champ coercitive (20°C)
mua = Br0/(mu0*Hcb0);                                                      % Permeabilité relative de l'AP
rho_cu = 1.72*1e-8;                                                        % Résisitivité du cuivre (ohm.m)
mvfer = 7800;                                                              % masse volumique du fer (kg/m3)
q = 1;                                                                     % Pertes fers spécifiques (W/kg)

%%  materiaux

mi_addmaterial('cu',1,1,0,0,0,0, 0, 0, 0, 0, 0,1,0);
mi_addmaterial('air' ,1,1,0,0,0,0,0,1,0,0,0);
mi_addmaterial('Fer lineaire' ,1000,1000,0,0,0,0,0,1,0,0,0);
mi_addmaterial('NdFeB' ,mua,mua,Hcb0,0,0,0,0,1,0,0,0);

%%  circuits
mi_addcircprop('a+',0,1);
mi_addcircprop('a-',0,1);
mi_addcircprop('b+',0,1);
mi_addcircprop('b-',0,1);
mi_addcircprop('c+',0,1);
mi_addcircprop('c-',0,1);

%% Conditions aux limites
mi_addboundprop('zero', 0, 0, 0, 0, 0, 0, 0, 0, 0);
mi_addboundprop('ps1', 0, 0, 0, 0, 0, 0, 0, 0, 4);
mi_addboundprop('pr1', 0, 0, 0, 0, 0, 0, 0, 0, 4);
mi_addboundprop('pr2', 0, 0, 0, 0, 0, 0, 0, 0, 4);
mi_addboundprop('pe1', 0, 0, 0, 0, 0, 0, 0, 0, 4);
mi_addboundprop('pe2', 0, 0, 0, 0, 0, 0, 0, 0, 4);
mi_addboundprop('pe3', 0, 0, 0, 0, 0, 0, 0, 0, 4);
mi_addboundprop('pe4', 0, 0, 0, 0, 0, 0, 0, 0, 4);
mi_addboundprop('pe5', 0, 0, 0, 0, 0, 0, 0, 0, 4);
mi_addboundprop('pe6', 0, 0, 0, 0, 0, 0, 0, 0, 4);

%% *** definition de la Géométrie

%% Linearization - 2D <- 3D

Req = (Rext+Rint)/2;                                                       % Rayon équivalent
lm = 2*pi*Req;                                                             % Longuer équivalente
prof_eq = (Rext^2-Rint^2)/(2*Req);                                         % Profondeur équivalente

deltay = 0;                                                                % Profondeur des aimants par rapport à la superfice du rotor
deltax = (lm/(2*P) - 2*pi*Req*(angle_aimant/360));                         % Distance entre deux aimants consécutifs
%deltax = 2;  %[mm]                                                        % Distance entre deux aimants consécutifs
delta_aimant = (2*pi*Req*(angle_aimant/360)+deltax);                       % Distance entre deux points équivalents de deux aimants consécutifs
delta_dents =  (volume_dents*lm/(3*P));                                    % longueur d'un dent
delta_encoche = (volume_encoche*lm/(3*P));                                 % longueur d'un encoche

%% Rotor

NodeR(1,:) = [0 h_e+h_S/2];                                                % point en bas à gauche
NodeR(2,:) = [0 h_e+h_S/2+h_R];                                            % point en haut à gauche
NodeR(3,:) = [lm h_e+h_S/2];                                               % point en bas à droite
NodeR(4,:) = [lm h_e+h_S/2+h_R];                                           % point en haut à droite

mi_drawrectangle(NodeR(1,:),NodeR(4,:));
mi_selectsegment((NodeR(1,:)+NodeR(2,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,1) ; mi_clearselected();
mi_selectsegment((NodeR(1,:)+NodeR(3,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,1) ; mi_clearselected();
mi_selectsegment((NodeR(2,:)+NodeR(4,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,1) ; mi_clearselected();
mi_selectsegment((NodeR(3,:)+NodeR(4,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,1) ; mi_clearselected();

%  Volume Rotor

hauteur_rotor = (1e-3)*(NodeR(2,2)-NodeR(1,2));
largeur_rotor = (1e-3)*(NodeR(3,1)-NodeR(1,1));
Vol_Rotor = (1e-3)*prof_eq*hauteur_rotor*largeur_rotor;

%--------------------------------------------------------------------------

% Basis of magnet

NodeAIM(1,:) = [deltax/2 h_e+h_S/2+deltay];                                              % point en bas à gauche
NodeAIM(2,:) = [deltax/2 h_e+h_S/2+deltay+l_aimant];                                     % point en haut à gauche
NodeAIM(3,:) = [deltax/2+2*pi*Req*(angle_aimant/360) h_e+h_S/2+deltay];                  % point en bas à droite
NodeAIM(4,:) = [deltax/2+2*pi*Req*(angle_aimant/360) h_e+h_S/2+deltay+l_aimant];         % point en haut à droite

% Drawing of magnet

mi_drawrectangle(NodeAIM(1,:),NodeAIM(4,:));
mi_selectsegment((NodeAIM(1,:)+NodeAIM(2,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,2) ; mi_clearselected();
mi_selectsegment((NodeAIM(1,:)+NodeAIM(3,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,2) ; mi_clearselected();
mi_selectsegment((NodeAIM(2,:)+NodeAIM(4,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,2) ; mi_clearselected();
mi_selectsegment((NodeAIM(3,:)+NodeAIM(4,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,2) ; mi_clearselected();

mi_selectgroup(2);
mi_copytranslate2(delta_aimant,0,2*P-1,4);

mi_selectgroup(1);
mi_selectgroup(2);
mi_mirror(0,0,lm,0)

%  Volume Aimants

hauteur_aimants = (1e-3)*(NodeAIM(2,2)-NodeAIM(1,2));
largeur_aimants = (1e-3)*(NodeAIM(3,1)-NodeAIM(1,1));
Vol_aimants = (1e-3)*2*P*prof_eq*hauteur_aimants*largeur_aimants;

% Masse Fer Culasse Rotorique
Mfer = mvfer*(Vol_Rotor-Vol_aimants);

%% Stator

% Bobines

NodeSB(1,:) = [-0 -h_b/2];                                                 % point en bas à gauche
NodeSB(2,:) = [-0  h_b/2];                                                 % point en haut à gauche
NodeSB(3,:) = [delta_encoche*0.5  -h_b/2];                                 % point en bas à droite
NodeSB(4,:) = [delta_encoche*0.5  h_b/2];                                  % point en haut à droite

mi_drawrectangle(NodeSB(1,:),NodeSB(4,:));
mi_selectsegment((NodeSB(1,:)+NodeSB(2,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,3) ; mi_clearselected();
mi_selectsegment((NodeSB(1,:)+NodeSB(3,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,3) ; mi_clearselected();
mi_selectsegment((NodeSB(2,:)+NodeSB(4,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,3) ; mi_clearselected();
mi_selectsegment((NodeSB(3,:)+NodeSB(4,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,3) ; mi_clearselected();

%Dents

NodeSD(1,:) = [delta_encoche*0.5 -h_S/2];                                  % point en bas à gauche
NodeSD(2,:) = [delta_encoche*0.5  h_S/2];                                  % point en haut à gauche
NodeSD(3,:) = [delta_encoche*0.5+delta_dents -h_S/2];                      % point en bas à droite
NodeSD(4,:) = [delta_encoche*0.5+delta_dents  h_S/2];                      % point en haut à droite

mi_drawrectangle(NodeSD(1,:),NodeSD(4,:));
mi_selectsegment((NodeSD(1,:)+NodeSD(2,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,4) ; mi_clearselected();
mi_selectsegment((NodeSD(1,:)+NodeSD(3,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,4) ; mi_clearselected();
mi_selectsegment((NodeSD(2,:)+NodeSD(4,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,4) ; mi_clearselected();
mi_selectsegment((NodeSD(3,:)+NodeSD(4,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,4) ; mi_clearselected();
mi_selectsegment((NodeSB(4,:)+NodeSD(2,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,4) ; mi_clearselected();
mi_selectsegment((NodeSB(3,:)+NodeSD(1,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,4) ; mi_clearselected();

% Mirroir bobine

mi_selectgroup(3);
mi_mirror(0.5*(NodeSD(2,:)+NodeSD(4,:)),0.5*(NodeSD(1,:)+NodeSD(3,:)))

% Translate

mi_selectgroup(3);
mi_selectgroup(4);
mi_copytranslate2(delta_encoche+delta_dents,0,3*P-1,4);


%% Déplacement Rotor et Entrefer

mi_selectgroup(1);
mi_selectgroup(2);
pas_init=(lm/(pas*P)-38*pi*Req/180);
mi_movetranslate2(pas_init,0,4)
mi_clearselected()

NodeE(1,:) = [0 (h_b/2)];
NodeE(2,:) = [0 0.95*h_e+h_S/2];
NodeE(3,:) = [pas_init 0.95*h_e+h_S/2];
NodeE(4,:) = [pas_init h_e+h_S/2];

mi_addnode(NodeE(1,:))
mi_addnode(NodeE(2,:))
mi_addnode(NodeE(3,:))
mi_addnode(NodeE(4,:))

mi_addsegment(NodeE(1,:),NodeE(2,:))
mi_addsegment(NodeE(2,:),NodeE(3,:))
mi_addsegment(NodeE(3,:),NodeE(4,:))

NodeE(5,:) = [lm (h_b/2)];
NodeE(6,:) = [lm 0.95*h_e+h_S/2];
NodeE(7,:) = [lm+pas_init 0.95*h_e+h_S/2];
NodeE(8,:) = [lm+pas_init h_e+h_S/2];

mi_selectsegment((NodeE(1,:)+NodeE(2,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,5) ; mi_clearselected();
mi_selectsegment((NodeE(2,:)+NodeE(3,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,5) ; mi_clearselected();
mi_selectsegment((NodeE(3,:)+NodeE(4,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,5) ; mi_clearselected();

mi_addnode(NodeE(5,:))
mi_addnode(NodeE(6,:))
mi_addnode(NodeE(7,:))
mi_addnode(NodeE(8,:))

mi_addsegment(NodeE(5,:),NodeE(6,:))
mi_addsegment(NodeE(6,:),NodeE(7,:))
mi_addsegment(NodeE(7,:),NodeE(8,:))

mi_selectsegment((NodeE(5,:)+NodeE(6,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,5) ; mi_clearselected();
mi_selectsegment((NodeE(6,:)+NodeE(7,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,5) ; mi_clearselected();
mi_selectsegment((NodeE(7,:)+NodeE(8,:))*0.5); mi_setsegmentprop('<None>', 0, 1, 0,5) ; mi_clearselected();

mi_selectgroup(5);
mi_mirror(0,0,lm,0);

%% *** Affectation des matériaux

% ======================================================
% Rotor

cord = mean(NodeR);
mi_addblocklabel(cord(1)+pas_init,cord(2));
mi_selectlabel(cord(1)+pas_init,cord(2));
mi_setblockprop('Fer lineaire',1,0,'Fer lineaire',0,1,0);
mi_clearselected( );
mi_selectlabel(cord(1)+pas_init,cord(2));
mi_copytranslate(0,-(2*(h_e+deltay)+h_S+h_R),1)
mi_clearselected();

% ======================================================
% Bloque aimant

direction = 1;
cord = mean(NodeAIM);
for i=1:2*P
    mi_addblocklabel(cord(1)+pas_init,cord(2));
    mi_selectlabel(cord(1)+pas_init,cord(2));
    mi_setblockprop('NdFeB',1,0,'<None>',90*direction,2,0);
    mi_clearselected( );
    mi_selectlabel(cord(1),cord(2));
    mi_movetranslate(delta_aimant*(2*P-i),0);
    mi_clearselected();
    mi_selectlabel(cord(1)+pas_init+delta_aimant*(2*P-i),cord(2));
    mi_copytranslate(0,-(2*(h_e+deltay)+h_S+l_aimant),1)
    mi_clearselected();
    direction = direction*(-1);
end
clear direction

% ======================================================
% Bobines

%Bobine A+

cord = mean(NodeSB);
mi_addblocklabel(cord(1),cord(2));
mi_selectlabel(cord(1),cord(2));
mi_setblockprop('cu',0,2,'a+',0,3,nb_spires)
mi_clearselected( );
mi_selectlabel(cord(1),cord(2));
mi_copytranslate(3*delta_dents+3*delta_encoche,0,(P-1))
mi_clearselected();

%Bobine A-

cord = mean(NodeSB)+[delta_dents+0.5*delta_encoche,0];
mi_addblocklabel(cord(1),cord(2));
mi_selectlabel(cord(1),cord(2));
mi_setblockprop('cu',0,2,'a-',0,3,nb_spires)
mi_clearselected( );
mi_selectlabel(cord(1),cord(2));
mi_copytranslate(3*delta_dents+3*delta_encoche,0,(P-1))
mi_clearselected();

%Bobine B+

cord = mean(NodeSB)+[delta_dents+1*delta_encoche,0];
mi_addblocklabel(cord(1),cord(2));
mi_selectlabel(cord(1),cord(2));
mi_setblockprop('cu',0,2,'b+',0,3,nb_spires)
mi_clearselected( );
mi_selectlabel(cord(1),cord(2));
mi_copytranslate(3*delta_dents+3*delta_encoche,0,(P-1))
mi_clearselected();

%Bobine B-

cord = mean(NodeSB)+[2*delta_dents+1.5*delta_encoche,0];
mi_addblocklabel(cord(1),cord(2));
mi_selectlabel(cord(1),cord(2));
mi_setblockprop('cu',0,2,'b-',0,3,nb_spires)
mi_clearselected( );
mi_selectlabel(cord(1),cord(2));
mi_copytranslate(3*delta_dents+3*delta_encoche,0,(P-1))
mi_clearselected();

%Bobine C+

cord = mean(NodeSB)+[2*delta_dents+2*delta_encoche,0];
mi_addblocklabel(cord(1),cord(2));
mi_selectlabel(cord(1),cord(2));
mi_setblockprop('cu',0,2,'c+',0,3,nb_spires)
mi_clearselected( );
mi_selectlabel(cord(1),cord(2));
mi_copytranslate(3*delta_dents+3*delta_encoche,0,(P-1))
mi_clearselected();

%Bobine C-

cord = mean(NodeSB)+[3*delta_dents+2.5*delta_encoche,0];
mi_addblocklabel(cord(1),cord(2));
mi_selectlabel(cord(1),cord(2));
mi_setblockprop('cu',0,2,'c-',0,3,nb_spires)
mi_clearselected();
mi_selectlabel(cord(1),cord(2));
mi_copytranslate(3*delta_dents+3*delta_encoche,0,(P-1))
mi_clearselected();

% ======================================================
% Dents

cord = mean(NodeSD);
mi_addblocklabel(cord(1),cord(2));
mi_selectlabel(cord(1),cord(2));
mi_setblockprop('air',1,0,'<None>',0,4,0);
mi_clearselected( );
mi_selectlabel(cord(1),cord(2));
mi_copytranslate(delta_dents+delta_encoche,0,3*P-1)
mi_clearselected();

% ======================================================
% Entrefer

cord = mean(NodeE);
mi_addblocklabel(cord(1),cord(2));
mi_selectlabel(cord(1),cord(2));
mi_setblockprop('air',0,2,'<None>',0,5,0);mi_clearselected( );
mi_selectlabel(cord(1),cord(2));
mi_clearselected();

mi_addblocklabel(cord(1),-cord(2));
mi_selectlabel(cord(1),-cord(2));
mi_setblockprop('air',0,2,'<None>',0,5,0);mi_clearselected( );
mi_selectlabel(cord(1),cord(2));
mi_clearselected();

%--------------------------------------------------------------------------
%% *** Conditions limites

% Rotor superfices externes

mi_selectsegment((NodeR(2,:)+NodeR(4,:))*0.5+[pas_init,0]);
mi_selectsegment((NodeR(2,1)+NodeR(4,1))*0.5,-2*(h_R+h_e)-h_S+(NodeR(2,2)));
mi_setsegmentprop('zero',0,1,0,1)
mi_clearselected()

% Rotor cotés

mi_selectsegment((NodeR(1,:)+NodeR(2,:))*0.5+[pas_init,0]);
mi_selectsegment((NodeR(3,:)+NodeR(4,:))*0.5+[pas_init,0]);
mi_setsegmentprop('pr1',0,1,0,1)
mi_clearselected()

mi_selectsegment((NodeR(1,1)+NodeR(1,1))*0.5 + pas_init,-2*h_e-h_R-h_S+(NodeR(1,2)+NodeR(1,2))*0.5);
mi_selectsegment((NodeR(3,1)+NodeR(4,1))*0.5 + pas_init,-2*h_e-h_R-h_S+(NodeR(3,2)+NodeR(4,2))*0.5);
mi_setsegmentprop('pr2',0,1,0,1)
mi_clearselected()

% Stator cotés

mi_selectsegment((NodeSB(1,:)+NodeSB(2,:))*0.5);
mi_selectsegment((NodeSB(1,1)+NodeSB(2,1))*0.5+lm,(NodeSB(1,2)+NodeSB(2,2))*0.5);
mi_setsegmentprop('ps1',0,1,0,4)
mi_clearselected()

% Entrefer

mi_selectsegment((NodeE(1,:)+NodeE(2,:))*0.5);
mi_selectsegment((NodeE(5,:)+NodeE(6,:))*0.5);
mi_setsegmentprop('pe1',0,1,0,5);
mi_clearselected()

mi_selectsegment((NodeE(2,:)+NodeE(3,:))*0.5);
mi_selectsegment((NodeE(6,:)+NodeE(7,:))*0.5);
mi_setsegmentprop('pe2',0,1,0,6)
mi_clearselected()

mi_selectsegment((NodeE(3,:)+NodeE(4,:))*0.5);
mi_selectsegment((NodeE(7,:)+NodeE(8,:))*0.5);
mi_setsegmentprop('pe3',0,1,0,7)
mi_clearselected()

mi_selectsegment((NodeE(1,1)+NodeE(2,1))*0.5,-(NodeE(1,2)+NodeE(2,2))*0.5);
mi_selectsegment((NodeE(5,1)+NodeE(6,1))*0.5,-(NodeE(5,2)+NodeE(6,2))*0.5);
mi_setsegmentprop('pe4',0,1,0,8)
mi_clearselected()

mi_selectsegment((NodeE(2,1)+NodeE(3,1))*0.5,-(NodeE(2,2)+NodeE(3,2))*0.5);
mi_selectsegment((NodeE(6,1)+NodeE(7,1))*0.5,-(NodeE(6,2)+NodeE(7,2))*0.5);
mi_setsegmentprop('pe5',0,1,0,9)
mi_clearselected()

mi_selectsegment((NodeE(3,1)+NodeE(4,1))*0.5,-(NodeE(3,2)+NodeE(4,2))*0.5);
mi_selectsegment((NodeE(7,1)+NodeE(8,1))*0.5,-(NodeE(7,2)+NodeE(8,2))*0.5);
mi_setsegmentprop('pe6',0,1,0,10)
mi_clearselected()

%% Déplacement pour mettre en phase le courant et la tension à vide

%% Sauvegarde

%--------------------------------------------------------------------------
mi_saveas('MS_axiale_RSR_NS_bobinage_concentre_V10.fem');

% %--------------------------------------------------------------------------
mi_analyze()
mi_loadsolution();

%% grandeur

Times=Tperiod/pas;
dt=Times;

%--------------------------------------------------------------------------
dtta=2*pi/P; %radians angle mecanique correspndant a l'angle electrique de 2pi
pastranslation=lm/(P*pas);
t=0:Times:Tperiod;

%--------------------------------------------------------------------------
%% Calculation of FLux at no load conditions
%--------------------------------------------------------------------------

for k=1:pas
    
    ia(k) = 0;
    ib(k) = 0;
    ic(k) = 0;
    iA(k) = Ieff*sqrt(2)*sin(((2*pi).*f)*t(k));
    iB(k) = Ieff*sqrt(2)*sin(((2*pi).*f)*t(k)-2*pi/3);
    iC(k) = Ieff*sqrt(2)*sin(((2*pi).*f)*t(k)-4*pi/3);
    
    % pour rajouter le courant dans les differents circuits a chaque pas dentaire
    
    mi_modifycircprop('a+',1,ia(k));
    mi_modifycircprop('a-',1,-ia(k));
    mi_modifycircprop('b+',1,ib(k));
    mi_modifycircprop('b-',1,-ib(k));
    mi_modifycircprop('c+',1,ic(k));
    mi_modifycircprop('c-',1,-ic(k));
    %--------------------------------------------------------------------------
    mi_analyze();mi_loadsolution();mi_zoomnatural();
    %%-------------------------------------------------------------------------
    a1=mo_getcircuitproperties('a+');
    a2=mo_getcircuitproperties('a-');
    %--------------------------------------------------------------------------
    b1=mo_getcircuitproperties('b+');
    b2=mo_getcircuitproperties('b-');
    %--------------------------------------------------------------------------
    c1=mo_getcircuitproperties('c+');
    c2=mo_getcircuitproperties('c-');
    %--------------------------------------------------------------------------
    Fluxboba0(k)=a1(3)-a2(3);
    Fluxbobb0(k)=b1(3)-b2(3);
    Fluxbobc0(k)=c1(3)-c2(3);
    %--------------------------------------------------------------------------
    
    % Couple methode tenseurmaxwell
    mo_groupselectblock(3);
    mo_groupselectblock(4);
    Couple_Maxwell_vide(k)=mo_blockintegral(18);
    mo_clearblock();
    %--------------------------------------------------------------------------
    %rotation
    mi_selectgroup(1);
    mi_selectgroup(2);
    mi_selectgroup(7);
    mi_selectgroup(10);
    mi_movetranslate2(pastranslation,0,4)
    mi_clearselected()
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------
figure('name','Flux à vide');
plot(1:360/(pas):360,Fluxboba0); hold on
plot(1:360/(pas):360,Fluxbobb0,'r'); hold on
plot(1:360/(pas):360,Fluxbobc0,'g'); hold on
xlabel('Angle(°)');
xlabel('Angle(°)');
ylabel('Flux à vide(wb)');

figure; plot(1:360/(pas):360,Couple_Maxwell_vide,'r');hold on;grid on
MeanTorque_vide=mean(Couple_Maxwell_vide);


% Calcul Fem à vide
dx = pastranslation ;
Fem = ((diff(Fluxboba0)/dx)*lm)*(N*2*pi/60); % Voir si il a le même résultat que les autres.
FemInitA=diff(Fluxboba0)/dt;
FemInitB=diff(Fluxbobb0)/dt;
FemInitC=diff(Fluxbobc0)/dt;
%--------------------------------------------------------------------------

figure('name','Tension à vide et Courant');
yyaxis left
plot(1:360/(pas-1):360,FemInitA,'b');hold on;grid on;
plot(1:360/(pas-1):360,FemInitB,'r');hold on;grid on;
plot(1:360/(pas-1):360,FemInitC,'g');hold on;grid on;
ylabel('FEM (V) à vide');
yyaxis right
plot(1:360/(pas-1):360,iA(1:pas-1),'b*');hold on;grid on;
plot(1:360/(pas-1):360,iB(1:pas-1),'r*');hold on;grid on;
plot(1:360/(pas-1):360,iC(1:pas-1),'g*');hold on;grid on;
ylabel('Courant (A)');
xlabel('Angle(°)');

%--------------------------------------------------------------------------
legend('PhaseA','PhaseB','PhaseC')
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% Calculation of FLux at load conditions - Alimentation Bobines
%--------------------------------------------------------------------------

% Mettre le Rotor dans la position iniciale.
mi_selectgroup(1);
mi_selectgroup(2);
mi_selectgroup(7);
mi_selectgroup(10);
mi_movetranslate2(-k*pastranslation,0,4)
mi_clearselected()

for k=1:pas
    
    ia(k) = Ieff*sqrt(2)*sin(((2*pi).*f)*t(k));
    ib(k) = Ieff*sqrt(2)*sin(((2*pi).*f)*t(k)-2*pi/3);
    ic(k) = Ieff*sqrt(2)*sin(((2*pi).*f)*t(k)-4*pi/3);
    
    % pour rajouter le courant dans les differents circuits a chaque pas dentaire
    
    mi_modifycircprop('a+',1,ia(k));
    mi_modifycircprop('a-',1,-ia(k));
    mi_modifycircprop('b+',1,ib(k));
    mi_modifycircprop('b-',1,-ib(k));
    mi_modifycircprop('c+',1,ic(k));
    mi_modifycircprop('c-',1,-ic(k));
    %--------------------------------------------------------------------------
    mi_analyze();mi_loadsolution();mi_zoomnatural();
    %%-------------------------------------------------------------------------
    a1=mo_getcircuitproperties('a+');
    a2=mo_getcircuitproperties('a-');
    %--------------------------------------------------------------------------
    b1=mo_getcircuitproperties('b+');
    b2=mo_getcircuitproperties('b-');
    %--------------------------------------------------------------------------
    c1=mo_getcircuitproperties('c+');
    c2=mo_getcircuitproperties('c-');
    %--------------------------------------------------------------------------
    Fluxboba0(k)=a1(3)-a2(3);
    Fluxbobb0(k)=b1(3)-b2(3);
    Fluxbobc0(k)=c1(3)-c2(3);
    %--------------------------------------------------------------------------
    
    % Pertes Fer Rotor
    
    % Rotor 1
    Cord_R = mean(NodeR);
    milieur_dent = mo_getb(Cord_R(1)+pastranslation*k+pas_init,Cord_R(2)) ;
    Bx_R1(k) = milieur_dent(1) ;
    By_R1(k) = milieur_dent(2) ;
    B_R1(k) = ((Bx_R1(k).^2)+(By_R1(k).^2)) ;
    
    % Rotor 2
    milieur_dent = mo_getb(Cord_R(1)+pastranslation*k+pas_init,-Cord_R(2)) ;
    Bx_R2(k) = milieur_dent(1) ;
    By_R2(k) = milieur_dent(2) ;
    B_R2(k) = ((Bx_R2(k).^2)+(By_R2(k).^2)) ;
    
    % Couple methode tenseurmaxwell
    mo_groupselectblock(3);
    mo_groupselectblock(4);
    Couple_Maxwell_charge(k)=mo_blockintegral(18);
    mo_clearblock();
    %--------------------------------------------------------------------------
    %rotation
    mi_selectgroup(1);
    mi_selectgroup(2);
    mi_selectgroup(7);
    mi_selectgroup(10);
    mi_movetranslate2(pastranslation,0,4)
    mi_clearselected()
    %--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------
figure('name','Flux en charge');
plot(1:360/(pas):360,Fluxboba0); hold on
plot(1:360/(pas):360,Fluxbobb0,'r'); hold on
plot(1:360/(pas):360,Fluxbobc0,'g'); hold on
xlabel('Angle(°)');
xlabel('Angle(°)');
ylabel('Flux en charge(wb)');

figure; plot(1:360/(pas):360,Couple_Maxwell_charge,'r');hold on;grid on
MeanTorque_charge = mean(Couple_Maxwell_charge);


% Calcul Fem
dx = pastranslation ;
Fem = ((diff(Fluxboba0)/dx)*lm)*(N*2*pi/60);
FemInitA=diff(Fluxboba0)/dt; 
FemInitB=diff(Fluxbobb0)/dt;
FemInitC=diff(Fluxbobc0)/dt;
%--------------------------------------------------------------------------

figure('name','Tension et Courant');
yyaxis left
plot(1:360/(pas-1):360,FemInitA,'b');hold on;grid on;
plot(1:360/(pas-1):360,FemInitB,'r');hold on;grid on;
plot(1:360/(pas-1):360,FemInitC,'g');hold on;grid on;
ylabel('FEM (V)');
yyaxis right
plot(1:360/(pas-1):360,ia(1:pas-1),'b*');hold on;grid on;
plot(1:360/(pas-1):360,ib(1:pas-1),'r*');hold on;grid on;
plot(1:360/(pas-1):360,ic(1:pas-1),'g*');hold on;grid on;
ylabel('Courant (A)');
xlabel('Angle(°)');

%--------------------------------------------------------------------------
legend('PhaseA','PhaseB','PhaseC')
%--------------------------------------------------------------------------

%% Rendement/ Pertes

% Géometrie Dents et enconches
largeur_encoche = (1e-3)*abs(NodeSB(3,1)-NodeSB(1,1));
hauteur_encoche = (1e-3)*abs(NodeSB(1,2)-NodeSB(2,2));
Sencoche = hauteur_encoche*largeur_encoche;                                % Surface d'un encoche
largeur_dent = abs(NodeSD(3,1)-NodeSD(1,1));

lcuivre = 2*(1e-3)*(prof_eq + largeur_dent);
Sfil_cu = Sencoche*Coeff_remplissage;
Rphase = rho_cu*(nb_spires^2*lcuivre)/(Sfil_cu);                           % Résistance par bobine

% Pertes Joules
P_Joules = 3*P*Rphase*Ieff^2 ;                                             % Machine triphasé (x3) et P Bobines chaque phase pour P paires de pôles (xP)

% Pertes Fer
Pf_R1 = q*Mfer*(f/50)*max(B_R1);
Pf_R2 = q*Mfer*(f/50)*max(B_R2);

% Total
Pf_Total = Pf_R2 + Pf_R1;

% Pertes Totales
P_total = Pf_Total + P_Joules;

%--------------------------------------------------------------------------
%% Calculation of FLux at load conditions - Méthode Simplifié
%--------------------------------------------------------------------------

% Mettre le Rotor dans la position iniciale.
mi_selectgroup(1);
mi_selectgroup(2);
mi_selectgroup(7);
mi_selectgroup(10);
mi_movetranslate2(-k*pastranslation,0,4)
mi_clearselected()

% Méthode Simplifié
temps = [0 pi/12 pi/6 pi/4]/(P);                                           % Angles par rapport à la partie électrique
deplacement = temps.*lm;                                                   % Déplacement mécanique par rapport aux angles électriques

for t=1:length(temps)
    
    
    
    %rotation
    mi_selectgroup(1);
    mi_selectgroup(2);
    mi_selectgroup(7);
    mi_selectgroup(10);
    mi_movetranslate2(pastranslation,0,4)
    mi_clearselected()
    
    
    
end

%% *** End :))))
   ia(k) = 0;
    ib(k) = 0;
    ic(k) = 0;
    iA(k) = Ieff*sqrt(2)*sin(((2*pi).*f)*t(k));
    iB(k) = Ieff*sqrt(2)*sin(((2*pi).*f)*t(k)-2*pi/3);
    iC(k) = Ieff*sqrt(2)*sin(((2*pi).*f)*t(k)-4*pi/3);
    
    % pour rajouter le courant dans les differents circuits a chaque pas dentaire
    
    mi_modifycircprop('a+',1,ia(k));
    mi_modifycircprop('a-',1,-ia(k));
    mi_modifycircprop('b+',1,ib(k));
    mi_modifycircprop('b-',1,-ib(k));
    mi_modifycircprop('c+',1,ic(k));
    mi_modifycircprop('c-',1,-ic(k));
    %--------------------------------------------------------------------------
    mi_analyze();mi_loadsolution();mi_zoomnatural();
    %%-------------------------------------------------------------------------
    a1=mo_getcircuitproperties('a+');
    a2=mo_getcircuitproperties('a-');
    %--------------------------------------------------------------------------
    b1=mo_getcircuitproperties('b+');
    b2=mo_getcircuitproperties('b-');
    %--------------------------------------------------------------------------
    c1=mo_getcircuitproperties('c+');
    c2=mo_getcircuitproperties('c-');
    %--------------------------------------------------------------------------
    Fluxboba0(k)=a1(3)-a2(3);
    Fluxbobb0(k)=b1(3)-b2(3);
    Fluxbobc0(k)=c1(3)-c2(3);
    %--------------------------------------------------------------------------
    
    % Couple methode tenseurmaxwell
    mo_groupselectblock(3);
    mo_groupselectblock(4);
    Couple_Maxwell_vide(k)=mo_blockintegral(18);
    mo_clearblock();
    %--------------------------------------------------------------------------

%--------------------------------------------------------------------------





