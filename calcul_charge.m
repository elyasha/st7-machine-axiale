%% Calculation of FLux at load conditions - Alimentation Bobines
%--------------------------------------------------------------------------
opendocument('MS_axiale_RSR_NS_bobinage_concentre.fem');
mi_saveas('temp_charge.fem');
mi_analyze()
mi_loadsolution();

Times=Tperiod/pas;
dt=Times;

%--------------------------------------------------------------------------
dtta=2*pi/P; %radians angle mecanique correspndant a l'angle electrique de 2pi
pastranslation=lm/(P*pas);
t=0:Times:Tperiod;

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
plot(1:360/(pas):360,Fluxboba0); hold on;grid on;
plot(1:360/(pas):360,Fluxbobb0,'r'); hold on;grid on;
plot(1:360/(pas):360,Fluxbobc0,'g'); hold on;grid on;
xlim([0 360]);
xticks(0:60:360)
xlabel('Angle(°)');
ylabel('Flux en charge(wb)');
title ('Flux en charge');

figure('name','Couple en charge');
plot(1:360/(pas):360,Couple_Maxwell_charge,'r');hold on;grid on
xlim([0 360]);
xticks(0:60:360)
MeanTorque_charge = mean(Couple_Maxwell_charge);
xlabel('Angle(°)');
ylabel('Couple (N.mm)');
str = ["Couple Moyen en Charge: ", string(MeanTorque_charge),"N.m"];
title (join(str));

% Calcul Fem
dx = pastranslation*1e-3 ;
FemA = ((diff(Fluxboba0)/dx)*(lm*1e-3)*(N/60));
FemB = ((diff(Fluxbobb0)/dx)*(lm*1e-3)*(N/60));
FemC = ((diff(Fluxbobc0)/dx)*(lm*1e-3)*(N/60));
%--------------------------------------------------------------------------

figure('name','Tension et Courant');
plot(1:360/(pas-1):360,FemA,'b');hold on;grid on;
plot(1:360/(pas-1):360,FemB,'r');hold on;grid on;
plot(1:360/(pas-1):360,FemC,'g');hold on;grid on;
%--------------------------------------------------------------------------
legend('PhaseA','PhaseB','PhaseC')
%--------------------------------------------------------------------------

%% Rendement/ Pertes

% Pertes Joules
P_Joules = 3*P*Rphase*Ieff^2 ;                                             % Machine triphasé (x3) et P Bobines chaque phase pour P paires de pôles (xP)

% Pertes Fer
Pf_R1 = q*Mfer*(f/50)*max(B_R1);
Pf_R2 = q*Mfer*(f/50)*max(B_R2);

% Total
Pf_Total = Pf_R2 + Pf_R1;

% Pertes Totales
P_total = Pf_Total + P_Joules;

% Densité de courant
J = sqrt(2)*Ieff/(Sfil_cu/nb_spires);                                      % Densité de courant [A/m²]
%--------------------------------------------------------------------------

% Output
 Couple_charge = Couple_Maxwell_charge;
 Couple_moyen(2) = MeanTorque_charge;
 Pertes (1) = P_total;


