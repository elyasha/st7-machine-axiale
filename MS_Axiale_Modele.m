function [Couple_vide,Couple_charge, Couple_M2,Couple_moyen, Pertes, Masse,Femmax_vide, J] = MS_Axiale_Modele(Parameters, Variables, pas)

%% Open FEMM, Open New Document, Define the magnetostatic problem

openfemm();
newdocument(0);
mi_probdef(0,'millimeters','planar', 1E-8, 60, 30);

%% Initialisation des sorties
Couple_charge = 0 ;
Couple_vide = 0;
Couple_charge = 0;
Couple_M2 = 0;
Couple_moyen = 0;
Pertes = 0;
Masse = 0;
Femmax_vid = 0;
J = 0;

%% Variables

Rint      = Variables(1);
Req       = Variables(2);
lm        = Variables(3);
prof_eq   = Variables(4);
h_R       = Variables(5);
h_S       = Variables(6);
h_e       = Variables(7);
h_b       = Variables(8);
P         = Variables(9);
nb_spires = Variables(10);

%% Parameters

%L                 = Parameters(1);
Rext              = Parameters(2);
volume_encoche    = Parameters(3);
volume_dents      = Parameters(4);
Coeff_remplissage = Parameters(5);
%Veff              = Parameters(6);
Ieff              = Parameters(7);
%Imax_pique        = Parameters(8);
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

%-----------------------------------
%% Dessin la machine et la teste
%-----------------------------------
%% Définition des paramétres 

% Géometrie aimants 
Surface_Rotor = (Rext^2-Rint^2)*pi;
Surface_aimant = 2*P*(pi*R_aimant^2);
Surf_rot_aimant = Surface_aimant/Surface_Rotor;                            % Rapport de la surface du rotor occupé par les aimants.
angle_aimant=((2*pi*Surf_rot_aimant)/(2*P))*(180/pi);                      % Surf_rot_aimant*100 % de la surface du rotor couverte par des aimants.

deltay = 0;                                                                % Profondeur des aimants par rapport à la superfice du rotor
deltax = (lm/(2*P) - 2*pi*Req*(angle_aimant/360));                         % Distance entre deux aimants consécutifs
delta_aimant = (2*pi*Req*(angle_aimant/360)+deltax);                       % Distance entre deux points équivalents de deux aimants consécutifs
delta_dents =  (volume_dents*lm/(3*P));                                    % longueur d'un dent
delta_encoche = (volume_encoche*lm/(3*P));                                 % longueur d'un encoche

% Champ Tournante 
f=P*N/60;        % [Hz]                                                    % Frequence champ tournante
Tperiod=1/f;     % [s]                                                     % Periode champ tournante
w=2*pi*f;        % [rad/s]                                                 % Vitesse champ tournante


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

% Géometrie et Volume Dents et enconches
largeur_encoche = (1e-3)*abs(NodeSB(3,1)-NodeSB(1,1));
hauteur_encoche = (1e-3)*abs(NodeSB(1,2)-NodeSB(2,2));
Sencoche = hauteur_encoche*largeur_encoche;                                % Surface d'un encoche
largeur_dent = abs(NodeSD(3,1)-NodeSD(1,1));
Vol_Dents = (1e-9)*3*P*largeur_dent*h_S*prof_eq; 
Vol_Bobines = 6*P*(1e-3)*Sencoche*prof_eq;  

% Caracteristiques électriques bobines
lcuivre = 2*(1e-3)*(prof_eq + largeur_dent);
Sfil_cu = Sencoche*Coeff_remplissage;
Rphase = rho_cu*(nb_spires^2*lcuivre)/(Sfil_cu);                           % Résistance par bobine


%% Déplacement Rotor et Entrefer
mi_selectgroup(1);
mi_selectgroup(2);
Posit_init = mean(NodeSD) - mean(NodeAIM);
pas_init= Posit_init(1);
mi_movetranslate2(pas_init,0,4)
mi_clearselected()


NodeE(1,:) = [0 (h_b/2)];
NodeE(2,:) = [0 0.80*h_e+h_S/2];
NodeE(3,:) = [pas_init 0.80*h_e+h_S/2];
NodeE(4,:) = [pas_init h_e+h_S/2];

mi_addnode(NodeE(1,:))
mi_addnode(NodeE(2,:))
mi_addnode(NodeE(3,:))
mi_addnode(NodeE(4,:))

mi_addsegment(NodeE(1,:),NodeE(2,:))
mi_addsegment(NodeE(2,:),NodeE(3,:))
mi_addsegment(NodeE(3,:),NodeE(4,:))

NodeE(5,:) = [lm (h_b/2)];
NodeE(6,:) = [lm 0.80*h_e+h_S/2];
NodeE(7,:) = [lm+pas_init 0.80*h_e+h_S/2];
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


%% Calcul de la masse

Maimants = Vol_aimants*mvaim;
Mdents = Vol_Dents*mvplastique;
Mbobines = Vol_Bobines*mvcu*Coeff_remplissage;
Masse = Mfer + Maimants+Mdents + Mbobines;

%% Sauvegarde

%--------------------------------------------------------------------------
mi_saveas('MS_axiale_RSR_NS_bobinage_concentre.fem');

mi_close;
% %--------------------------------------------------------------------------

%Calcul_vide; 
calcul_charge_4pts;
end






