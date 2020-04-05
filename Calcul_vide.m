opendocument('MS_axiale_RSR_NS_bobinage_concentre.fem');
mi_saveas('temp_vide.fem');
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
    
    ia(k) = Ieff*sqrt(2)*sin(((2*pi).*f)*t(k));
    ib(k) = Ieff*sqrt(2)*sin(((2*pi).*f)*t(k)-2*pi/3);
    ic(k) = Ieff*sqrt(2)*sin(((2*pi).*f)*t(k)-4*pi/3);
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
mi_analyze()
mi_loadsolution();
end
%--------------------------------------------------------------------------
figure('name','Flux à vide');
plot(1:360/(pas):360,Fluxboba0); hold on;grid on;
plot(1:360/(pas):360,Fluxbobb0,'r'); hold on;grid on;
plot(1:360/(pas):360,Fluxbobc0,'g'); hold on;grid on;
xlim([0 360]);
xticks(0:60:360)
xlabel('Angle(°)');
ylabel('Flux à vide(wb)');
title ('Flux à vide');

figure('name','Couple à vide');
plot(1:360/(pas):360,Couple_Maxwell_vide,'r');hold on;grid on
xlim([0 360]);
xticks(0:60:360)
MeanTorque_vide = mean(Couple_Maxwell_vide);
xlabel('Angle(°)');
ylabel('Couple (N.mm)');
str = ["Couple Moyen en Charge: ", string(MeanTorque_charge),"N.m"];
title (join(str));

% Calcul Fem à vide
dx = pastranslation*1e-3;
FemA = ((diff(Fluxboba0)/dx)*(lm*1e-3)*(N/60));
FemB = ((diff(Fluxbobb0)/dx)*(lm*1e-3)*(N/60));
FemC = ((diff(Fluxbobc0)/dx)*(lm*1e-3)*(N/60));

figure('name','Tension et Courant');
yyaxis left
plot(1:360/(pas-1):360,ia,'b');hold on;grid on;
plot(1:360/(pas-1):360,ib,'r');hold on;grid on;
plot(1:360/(pas-1):360,ic,'g');hold on;grid on;
ylabel('FEM à vide (V)');
yyaxis right
plot(1:360/(pas-1):360,FemA,'b');hold on;grid on;
plot(1:360/(pas-1):360,FemB,'r');hold on;grid on;
plot(1:360/(pas-1):360,FemC,'g');hold on;grid on;
ylabel('Courant  (A)');
legend('PhaseA','PhaseB','PhaseC')
xlim([0 360]);
xticks(0:60:360)
xlabel('Angle(°)');
title ('FEM à vide et courant');
%--------------------------------------------------------------------------
 
% Output

Couple_vide = Couple_Maxwell_vide;

Femmax_vide = max(FemA);

mi_close;

%--------------------------------------------------------------------------
