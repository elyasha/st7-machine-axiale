%% Calculation of FLux at load conditions - Méthode Simplifié
%--------------------------------------------------------------------------
opendocument('MS_axiale_RSR_NS_bobinage_concentre.fem');
mi_saveas('temp_4ptcharge.fem');
mi_analyze()
mi_loadsolution();

Times=Tperiod/pas;
dt=Times;

%--------------------------------------------------------------------------
dtta=2*pi/P; %radians angle mecanique correspndant a l'angle électrique de 2pi
pastranslation=lm/(P*pas);
t=0:Times:Tperiod;

mi_analyze()
mi_loadsolution();

% Méthode Simplifiée

tempresolution = [0 pi/12 pi/6 pi/4]*Tperiod/(2*pi);
DegElec = tempresolution*2*pi/Tperiod/P;

deplacement_rotor = DegElec*Req;

for k=1:length(tempresolution)
    
    ia(k) = Ieff*sqrt(2)*sin(((2*pi).*f)*tempresolution(k));
    ib(k) = Ieff*sqrt(2)*sin(((2*pi).*f)*tempresolution(k)-2*pi/3);
    ic(k) = Ieff*sqrt(2)*sin(((2*pi).*f)*tempresolution(k)-4*pi/3);
    
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
    Bx_R1_M2(k) = milieur_dent(1) ;
    By_R1_M2(k) = milieur_dent(2) ;
    B_R1_M2(k) = ((Bx_R1_M2(k).^2)+(By_R1_M2(k).^2)) ;
    
    % Rotor 2
    milieur_dent = mo_getb(Cord_R(1)+pastranslation*k+pas_init,-Cord_R(2)) ;
    Bx_R2_M2(k) = milieur_dent(1) ;
    By_R2_M2(k) = milieur_dent(2) ;
    B_R2_M2(k) = ((Bx_R2_M2(k).^2)+(By_R2_M2(k).^2)) ;
    
    % Couple methode tenseurmaxwell
    mo_groupselectblock(3);
    mo_groupselectblock(4);
    Couple_Maxwell_charge_M2(k)=mo_blockintegral(18);
    mo_clearblock();
    %--------------------------------------------------------------------------
    %rotation
    mi_selectgroup(1);
    mi_selectgroup(2);
    mi_selectgroup(7);
    mi_selectgroup(10);
    mi_movetranslate2(deplacement_rotor(k),0,4);
    mi_clearselected()
    
end

MeanTorque_charge_M2 = mean(Couple_Maxwell_charge_M2);

% Rendement/ Pertes

% Pertes Joules
P_Joules_M2 = 3*P*Rphase*Ieff^2 ;                                          % Machine triphasé (x3) et P Bobines chaque phase pour P paires de pôles (xP)

% Pertes Fer
Pf_R1_M2 = q*Mfer*(f/50)*max(B_R1_M2);
Pf_R2_M2 = q*Mfer*(f/50)*max(B_R2_M2);

% Total
Pf_Total_M2 = Pf_R2_M2 + Pf_R1_M2;

% Pertes Totales
P_total_M2 = Pf_Total_M2 + P_Joules_M2;

%% Output

Couple_M2       = Couple_Maxwell_charge_M2; 
Couple_moyen = MeanTorque_charge_M2;
Pertes       = P_total_M2;
mi_close;
