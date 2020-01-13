clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Partie I  :  Etude Longitudinale   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Mesure Directe
Coef_Conversion = 50/100; %% cm à l'échelle 1/50ème vers m à l'échelle 1
Coef_DegRad = pi/180;
Coef_RadDeg = 57.3;

L_Ref = 16.3;

Envergure = 9;
Envergure_ext = 7;
Corde_emplanture = 13;
CMA = 8.7;

d_x0_CMA = 16.5;
d_x0_Corde_emplanture = 12;

Fleche_50_deg = 36
Fleche_0_deg = 60
Fleche_50_rad = Fleche_50_deg*Coef_DegRad;
Fleche_0_rad = Fleche_0_deg*Coef_DegRad;

L_Fu = 4;
h_Fu = 13;

Base_al = 13;
Hauteur_al = 7;

Rayon_Nez = 1.6;          %% Ici égales
Rayon_EntAir = 1.6;       %%     //


%% Dimention & Surface réelle (en m²)

L_Ref_reelle = L_Ref*Coef_Conversion

S_v0_Ref = (L_Ref*Envergure)*Coef_Conversion^2

S_v0_Fu = (L_Fu * h_Fu)*Coef_Conversion^2 
S_v0_Ext = (Base_al * Hauteur_al)*Coef_Conversion^2

S_nez = pi*((Rayon_Nez*Coef_Conversion)^2)
S_EntAir = pi*((Rayon_EntAir*Coef_Conversion)^2)

S_alaire = S_v0_Ref+S_nez+S_EntAir


%% Calcul de Cza_i (en nbr/rad)

    %% Donnée 
    Cza_EntAir = 2;
    Cza_Fu = 2;
    
    %% Calcul de l'allongement
    Lambda = ((Envergure_ext*Coef_Conversion)^2) / (Corde_emplanture*Coef_Conversion)
    
    %% Calcul Cza_voilure
    Cza_voilure = (2*pi*Lambda)/(2 + sqrt(4+((Lambda^2) / (cos(Fleche_50_rad)^2))))
    
    %% Calcul Cza_interaction_fuselage
    Cza_interaction_fuselage = Cza_voilure  %%Hypothèse donnée
    
%%  Calcul de Cza_0 (en nbr/rad)  

    Cza_0 = (S_nez*Cza_Fu + S_EntAir*Cza_EntAir + S_v0_Ext*Cza_voilure + S_v0_Fu*Cza_interaction_fuselage)/S_v0_Ref

%%  Calcul de X_foyer (m à l'échelle 1) 
    
    %%Foyer Intermediaire
    X_ref = (53/100)*L_Ref*Coef_Conversion+8.5*Coef_Conversion;
    X_FuAv = 4.1*Coef_Conversion;
    X_EntAir = 9.9*Coef_Conversion;
    X_Voilure = ((((25+(Fleche_0_deg/10))/100)* (CMA*Coef_Conversion)))+(d_x0_CMA*Coef_Conversion);
    X_Interaction_Fuselage = ((25+(Fleche_0_deg/3.6))/100)* (Corde_emplanture*Coef_Conversion)+(d_x0_Corde_emplanture*Coef_Conversion);
    %%Foyer%%
    X_foyer = (S_nez*Cza_Fu*X_FuAv + S_EntAir*Cza_EntAir*X_EntAir + S_v0_Ext*Cza_voilure*X_Voilure + S_v0_Fu*Cza_interaction_fuselage*X_Interaction_Fuselage)/(S_nez*Cza_Fu + S_EntAir*Cza_EntAir + S_v0_Ext*Cza_voilure + S_v0_Fu*Cza_interaction_fuselage);
    
    %%Ratio en %de L_Ref
    X_ref_ratio= X_ref/L_Ref*100
    X_FuAv_ratio= X_FuAv/L_Ref*100
    X_EntAir_ratio= X_EntAir/L_Ref*100
    X_Voilure_ratio= X_Voilure/L_Ref*100
    X_Interaction_Fuselage_ratio= X_Interaction_Fuselage/L_Ref*100
    X_foyer_ratio= X_foyer/L_Ref*100
    
%% Etude de stabilité 

    %% Cm_alpha_CG 
    
    Cm_alpha_CG = -((X_foyer-X_ref)*Cza_0)/L_Ref
    
    %% MS_vol_rectiligne 
    
    MS_vol_rectiligne = -((X_ref-X_foyer)/L_Ref)*100  %%Pourcentage de Statibilité

    %% Question 4 
    
        %% Pour 4% de Stabilité on Cherche -((X_ref - X_foyer)/L_Ref)*100 = +4 
        %% Pour 4% d' Instabilité on Cherche -((X_ref - X_foyer)/L_Ref)*100 = -4 
        %% Pour un avion neutre on Cherche (Xcg - X_foyer) /L_Ref = 0 , d'où Xcg = X_foyer 

        syms X_ref_4
        S_4 = solve(-((X_ref_4 - X_foyer)/L_Ref)*100 == 4 ,X_ref_4);
        sol_4 = vpa(S_4)
        
        syms X_ref_moins4
        S_moins4 = solve(-((X_ref_moins4 - X_foyer)/L_Ref)*100 == -4 ,X_ref_moins4);
        sol_moins4 = vpa(S_moins4)
        
    %% Calcul du Cz_dm & Cm_dm

        %Coefficient
        Profondeur_Gouverne = (1.5/7.4)*100 %Pourcentage de la corde
        Kcor = 0.53;
        Ec = 0.3;
        DKenv=1 ;
        
        %%Calcul de Cz_d
        Khyper = Kcor*(1+Ec)*DKenv
        Cz_d = Cza_voilure * Khyper
        
        %%Calcul de Cz_dm
        Cz_dm = Cz_d * (S_v0_Ext/S_v0_Ref)
       
        %%Ratio de X_foyer_gouverne en pourcentage de L_Ref
        X_fg = 21.5*Coef_Conversion
       
        %%Calcul de Cm_dm à 53%
        Cm_dm = -((X_fg-X_ref)*Cz_dm)/(L_Ref_reelle)
   
        
        
        
 %% 2 - Equilibrage de l'avion

        %Donnée
        m = 9000; %%kg
        nz = 1; %% Facteur de charge en g 
        Cz0 = 0; %%cf p26
        Cm0 = 0.01; %%cf p26
        g=9.81;
        
        %%Incidence de braquage de gouverne à nz=1g
        rho = 0.65; %% kg/m^3
        V=253; %% m/s
        M = 0.8;
        z = 20000; %% ft
        pente = 0; %% deg
        
        qbar = (rho*V^2)/2;
        
        %% Résolution Système en vol
        syms alpha dm
        
        E1 = (qbar*S_v0_Ref)*(Cz0+(Cza_0*alpha)+(Cz_dm*dm))==m*g*cosd(pente);
        E2 = qbar*S_v0_Ref*L_Ref_reelle*(Cm0+Cm_alpha_CG*alpha+Cm_dm*dm)==0;
        
        [alpha,dm] = solve(E1,E2,alpha,dm);
        
        alpha_deg = vpa(alpha)*Coef_RadDeg
        dm_deg = vpa(dm)*Coef_RadDeg
        
        alpha_rad = vpa(alpha);
        dm_rad = vpa(dm);
        
                    %% NB : Avion presque neutre -> 2 valeurs positives
        
        %%Vitesse de l'avion et Braquage de la gouverne d'équilibre en approche à nz=1g
        incidence_approche = 14*Coef_DegRad; %% deg
        pente_approche = -3*Coef_DegRad; %% deg
        z_approche = 1000; %% ft
        rho_approche = 1.19;
           
         %% Résolution Système en approche
         
        syms V_approche dm_approche
        
        E1_approche = (((rho_approche*V_approche^2)/2)*S_v0_Ref)*(Cz0+(Cza_0*incidence_approche)+(Cz_dm*dm_approche))==m*g*cosd(pente_approche);
        E2_approche = ((rho_approche*V_approche^2)/2)*S_v0_Ref*L_Ref_reelle*(Cm0+Cm_alpha_CG*incidence_approche+Cm_dm*dm_approche)==0;
                
        [V_approche,dm_approche] = solve(E1_approche,E2_approche,V_approche,dm_approche);
        
        V_approche = vpa(V_approche)
        V_approche_noeud = V_approche*1.94384
        dm_deg_approche = vpa(dm_approche)*Coef_RadDeg
        dm_rad_approche = vpa(dm_approche)
        
        qbar_app = (rho*V_approche(2,1)^2)/2;
        
        %% Influence de la position CG (et donc marge Statique)sur les deux équilibre précédent
        p = 1;
        
        for X_ref_var = sol_4:0.1:sol_moins4
           
            MS_vol_rectiligne_var(1,p) = -((X_ref_var-X_foyer)/L_Ref)*100;
            
            %%Calcul de Cm_dm
            Cm_dm_var(1,p) = -((X_fg-X_ref_var)*Cz_dm)/(L_Ref_reelle);
            
            %%Calcul de Cm_alpha_CG 
            Cm_alpha_CG_var(1,p) = -((X_foyer-X_ref_var)*Cza_0)/L_Ref;
            
            
            %% Résolution Système en vol                                      REFAIRE RESOLUTION THEORIQUE 
                        
            Sec_Term = (m*g*cosd(pente_approche))/(qbar*S_v0_Ref);                 %%2nd membre de la premiere equation du systeme
            Coef_Dm(1,p) =(-Cm0-Cm_dm_var(1,p))/Cm_alpha_CG_var(1,p);           %%Constante des equation de resolution

            dm_rad_var(1,p)=((Sec_Term-Cz0)/(Coef_Dm(1,p)*Cza_0+Cz_dm));   
            alpha_rad_var(1,p) = (Coef_Dm(1,p) * dm_rad_var(1,p));      
            
            dm_deg_var(1,p) = dm_rad_var(1,p)*Coef_RadDeg;                      %%Converrsion
            alpha_deg_var(1,p) = alpha_rad_var(1,p)*Coef_RadDeg;                %%Conversion
            
            %% Résolution Système en approche                                 REFAIRE RESOLUTION THEORIQUE 
                        
            Sec_Term = m*g*cosd(pente_approche);                              %%2nd membre de la premiere equation du systeme
            Coef_Dm(1,p) =(-Cm0-Cm_alpha_CG_var(1,p)/Cm_dm_var(1,p));         %%Constante des equation de resolution

            dm_rad_var_app(1,p)=((Sec_Term-Cz0)/(Coef_Dm(1,p)*Cza_0+Cz_dm))/10^5; 
            V_approche_var(1,p) = sqrt(Sec_Term/((S_v0_Ref*(Cz0+Cza_0*incidence_approche+Cz_dm*dm_rad_var_app(1,p))))*rho_approche);       
            
            dm_deg_var_app(1,p) = dm_rad_var(1,p)*Coef_RadDeg;                    %%Conversion
            V_approche_var_noeud(1,p) = V_approche_var(1,p)*1.94384;                 %%Conversion

            p = p+1;

        end

        % Courbes phase en Vol
        figure(1)        
        plot(sol_4:0.1:sol_moins4,dm_deg_var)
        grid on
        title('Evolution de \delta_m en dégré selon variation de centrage en Vol')
        xlabel('X_ref (en m à echelle 1)')
        ylabel('Angle de gouverne (deg °)')
        
        figure(2)        
        plot(sol_4:0.1:sol_moins4,alpha_deg_var)
        grid on
        title('Evolution de \alpha en dégré selon variation de centrage en Vol')
        xlabel('X_ref (en m à echelle 1)')
        ylabel('Angle incidence (deg °)')
        
        %% Courbes phase d'approche
        figure(3)        
        plot(sol_4:0.1:sol_moins4,dm_deg_var_app)
        grid on
        title('Evolution de \delta_m en dégré selon variation de centrage en Approche')
        xlabel('X_ref (en m à echelle 1)')
        ylabel('Angle de gouverne (deg °)')
        
        figure(4)        
        plot(sol_4:0.1:sol_moins4,V_approche_var_noeud)
        grid on
        title('Evolution de la Vitesse approche en noeuds selon variation de centrage en Approche')
        xlabel('X_ref (en m à echelle 1)')
        ylabel('Vitesse (noeuds)')
        
        
        
        
%% 3 - Petits mouvement autour de l'équilibre

       %% Matrice d'inertie
       R_yy = 2.5; %% m 
       B_inertie = m*R_yy^2;
       
       %% PHASE DE VOL
       
               %% Mise sous forme matricielle X_point = A X + B U  (PHASE DE VOL)
               Cmq = -0.2; %% rad
               Fp = 0.5*rho*V^2*S_v0_Ref*Cza_0; %% Newton

               Mq = ((qbar*S_v0_Ref*L_Ref_reelle)/B_inertie) * (L_Ref_reelle/V) * Cmq;
               M_alpha = ((qbar*S_v0_Ref*L_Ref_reelle)/B_inertie) * Cm_alpha_CG;
               M_dm = ((qbar*S_v0_Ref*L_Ref_reelle)/B_inertie) * Cm_dm;

               Z_alpha = -((qbar*S_v0_Ref)/(m*V))*Cza_0 %%+ ((Fp/m*V)*cos(alpha_rad));
               Z_dm = -((qbar*S_v0_Ref)/(m*V))*Cz_dm;

               A = [ Mq , M_alpha ; 1 , Z_alpha]
               A = double(A)
               B = [ M_dm ; Z_dm]

               %% Calcul d'ammortissement et periode propre (PHASE DE VOL)

               W_naturelle = sqrt(Mq*Z_alpha - M_alpha);
               Amortissement = (-(Mq + Z_alpha))/(2 * W_naturelle)          %% Valeur aberrante, check utlisation coef
               W_propre = W_naturelle * sqrt(1-Amortissement^2);
               T_propre = vpa((2*pi)/W_propre)

       
       %% PHASE D'APPROCHE
       
               %% Mise sous forme matricielle X_point = A X + B U  (PHASE D'APPROCHE)
               Cmq_app = -0.2; %% rad
               Fp_app = 0.5*rho*V_approche(2,1)^2*S_v0_Ref*Cza_0; %% Newton

               Mq_app = ((qbar*S_v0_Ref*L_Ref_reelle)/B_inertie) * (L_Ref_reelle/V_approche(2,1)) * Cmq_app;
               M_alpha_app = ((qbar*S_v0_Ref*L_Ref_reelle)/B_inertie) * Cm_alpha_CG;
               M_dm_app = ((qbar*S_v0_Ref*L_Ref_reelle)/B_inertie) * Cm_dm;

               Z_alpha_app = -(((qbar*S_v0_Ref)/(m*V_approche(2,1)))*Cza_0) %%+ ((Fp_app/m*V_approche(2,1))*cos(incidence_approche));
               Z_dm_app = -((qbar*S_v0_Ref)/(m*V_approche(2,1)))*Cz_dm;

               A_app = [ Mq_app , M_alpha_app ; 1 , Z_alpha_app]
               A_app = double(A_app)
               B_app = [ M_dm_app ; Z_dm_app]
               B_app = double(B_app)


               %% Calcul d'ammortissement et periode propre (PHASE D'APPROCHE)

               W_naturelle_app = sqrt(Mq_app*Z_alpha_app - M_alpha_app);
               Amortissement_app = (-(Mq_app + Z_alpha_app))/(2 * W_naturelle_app)  %% Valeur Abérhente, check utlisation coef
               W_propre_app = W_naturelle_app * sqrt(1-Amortissement_app^2);
               T_propre_app = vpa((2*pi)/W_propre_app)

               
       %% Equation de sortie et calcul du facteur de charge au CG
       
               %% PHASE EN VOL
               
                   %% Calcul facteur de charge
                   nz = (qbar*S_v0_Ref*(Cza_0*alpha_rad+Cz_dm*dm_rad))/(m*g)

                        %% Y = C X + D U

                   %% Matrice C

                   C = [0 , 1 ; 1 , 0 ; 0 , qbar*S_v0_Ref*Cza_0/(m*g)]

                   %% Matrice D

                   D = [0 ;0 ;qbar*S_v0_Ref*Cz_dm/(m*g)]
                   
                   
                %% PHASE EN APPROCHE
               
                   %% Calcul facteur de charge
                   nz_app = (qbar_app*S_v0_Ref*(Cza_0*incidence_approche+Cz_dm*dm_rad_approche(1,1)))/(m*g)

                        %% Y = C X + D U

                   %% Matrice C

                   C_app = [0 , 1 ; 1 , 0 ; 0 , qbar_app*S_v0_Ref*Cza_0/(m*g)]
                   C_app = double(C_app)
                   %% Matrice D

                   D_app = [0 ;0 ;qbar_app*S_v0_Ref*Cz_dm/(m*g)]
                   D_app = double(D_app)
                   
        %% Comparaison
                    
               SYS_Vol=ss(A,B,C,D);
               damp(SYS_Vol)
               Amortissement
               SYS_App=ss(A_app,B_app,C_app,D_app);
               damp(SYS_App)
               Amortissement_app
               
                        %% On retrouve les même amortissement et période entre fonction Damp et Calcul Théorique a peu de chose près 
       
       %% Verification du calcul modal

                %% Plusieurs centrage 
                %% Calculer sur courbe Amortissement et T_propre
                
                sim('longi')
       
%                     %% Courbes phase en Vol
%                     
%                     figure()        
%                     plot(Temps,Commande_pilote)
%                     grid on
%                     title('Courbe Commande Pilote en Vol')
%                     xlabel('Temps (s)')
%                     ylabel('Commande')
%                     
%                     figure()        
%                     plot(Temps,profondeur_realisee)
%                     grid on
%                     title('Courbe de Profondeur \delta_m realisee en Vol')
%                     xlabel('Temps (s)')
%                     ylabel('Profondeur (°)')
%                     
%                     figure()        
%                     plot(Temps,incidence)
%                     grid on
%                     title('Courbe Incidence \alpha en Vol')
%                     xlabel('Temps (s)')
%                     ylabel('Incidence (°)')
%                     
%                     figure()        
%                     plot(Temps,vitesse_tangage)
%                     grid on
%                     title('Courbe de Vitesse de tangage en Vol')
%                     xlabel('Temps (s)')
%                     ylabel('Vitesse')
%                     
%                     figure()        
%                     plot(Temps,facteur_charge)
%                     grid on
%                     title('Courbe Facteur de charge Pilote en Vol')
%                     xlabel('Temps (s)')
%                     ylabel('N_z')                   %% BOUCLER
%                     
                 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Partie II  :  Etude Transversale   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calcul des coefficients transversaux


    
    %% Calcul de C_yB
        Lambda_derive =2;
        Fleche_derive_50_rad =38*Coef_DegRad;
        S_derive = 4.455;
        
        C_yb_pointeavant = 2 %%rad.s-1
        C_yb_entreeair = 2 %%rad.s-1
        C_yb_derive = (2*pi*Lambda_derive)/(2 + sqrt(4+((Lambda_derive^2) / (cos(Fleche_derive_50_rad)^2)))); %%rad.s-1
        
        Somme_C_yb = (S_nez*C_yb_pointeavant) + (S_EntAir*C_yb_entreeair) + (S_derive*C_yb_derive);
        
        C_yB = Somme_C_yb/S_v0_Ref
   
    %%%%%%%%%%%%%%%%%
    %   A REVOIR   %%
    %%%%%%%%%%%%%%%%%
    
    %% Calcul de C_nB

        x_PA = (2-X_ref)
        x_EntAir = (X_EntAir-X_ref)
        x_Derive = 9.7*Coef_Conversion
        
        Somme_C_yb_Xref = (S_nez*C_yb_pointeavant*x_PA) + (S_EntAir*C_yb_entreeair*x_EntAir) + (S_derive*C_yb_derive*x_Derive);
        
        C_nB = Somme_C_yb_Xref/(S_v0_Ref*L_Ref_reelle)
   
    %% Calcul de C_nr
        
        Cnr_d = ((C_yb_derive*x_Derive^2)/L_Ref_reelle^2)*(S_derive/S_v0_Ref)
        Cnr = 2*Cnr_d
    
    
    %% Coefficient transversaux
    
        C_yB = -0.45;
        C_nB = 0.15;
        
        C_lB_derive = -0.055;
        C_lB_fleche = -0.15;
        C_lB_diedre = -0.445;
        C_lB_hauteur_voilure = 0.05;
        
        C_lB = C_lB_derive + C_lB_fleche*(14/57.3)+C_lB_hauteur_voilure
        
        C_y_dn = 0.11
        C_l_dn = 0.025
        C_n_dn = -0.08
             
        C_y_dl = -0.15
        C_l_dl = 0.112
        C_n_dl = -0.15
        
        C_lp = -0.15
        C_nr = -0.23
        
        
%% Equilibraga avion

    
    %% Calcul du dérapage 
    V_vent = 15; %% Kts
    
    Beta_rad = asin(V_vent/V_approche_noeud(2,:)); %%rad
    Beta_deg = Beta_rad*Coef_RadDeg %%deg
    
    %% Resolution Système 
    second_terme_sys = ((m*g)/(qbar_app*S_v0_Ref))*cos(14/57.3);
    
    A_sys = [-second_terme_sys , C_y_dl  ,C_y_dn ; 0, C_l_dl , C_l_dn ; 0 ,C_n_dl , C_n_dn];
    
    
    B_sys = [-C_yB ; -C_lB ; -C_nB]*Beta_rad;
    
    X_sys = inv(A_sys)*B_sys;
    
    Phi = asin(X_sys(1,1));
    Phi_deg = Phi*Coef_RadDeg
    
    D_l0 = X_sys(2,1);
    D_l0_deg = D_l0*Coef_RadDeg
    
    D_n0 = X_sys(3,1); 
    D_n0_deg = D_n0*Coef_RadDeg
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Partie III  :  Systeme de Commande de Vol     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ms_ref = 0.1;

X_ref_ar = sol_moins4;
X_ref_av = sol_4;

    %%Calcul de G_alpha

    Cm_alpha_arriere = -((X_foyer-X_ref_ar)*Cza_0)/L_Ref;
    Cm_dm_arriere = -((X_fg-X_ref_ar)*Cz_dm)/(L_Ref_reelle);

    G_alpha=-((0.1*Cza_0+Cm_alpha_arriere)/(Cm_dm_arriere-0.1*Cz_dm))
    
    sim('longi.mdl')    
    
    [A_longi2,B_longi2,C_longi2,D_longi2]=linmod('longi2')
    [A_longi3,B_longi3,C_longi3,D_longi3]=linmod('longi3')
    
    damp(A_longi)
    Ms_sys = (B_inertie/(m*L_Ref_reelle*V))*(A_longi3(1,2)/A_longi3(2,2))

    %%Calcul de G_q

    Cm_alpha_avant = -((X_foyer-X_ref_av)*Cza_0)/L_Ref
    Cm_dm_avant = -((X_fg-X_ref_av)*Cz_dm)/(L_Ref_reelle)   
    
    M_alpha_avant=((qbar*S_v0_Ref*L_Ref_reelle)/B_inertie)*Cm_alpha_avant
    M_dm_avant=((qbar*S_v0_Ref*L_Ref_reelle)/B_inertie)*Cm_dm_avant
    
    syms G_q
    
    Z_a_chap=Z_alpha+Z_dm*G_alpha;
    Z_q_chap=1+Z_dm*G_q;
    M_a_chap=M_alpha_avant+M_dm_avant*G_alpha;
    M_q_chap=Mq+M_dm_avant*G_q;
         
    Q_q = solve((-(M_q_chap+Z_a_chap))/(2*sqrt(M_q_chap*Z_a_chap-M_a_chap*Z_q_chap)) == 0.7 ,G_q);
    Q_q_val = vpa(Q_q)
    
    [A_longi4,B_longi4,C_longi4,D_longi4]=linmod('longi4')
    
    damp(A_longi4)
    Ms_sys = (B_inertie/(m*L_Ref_reelle*V))*(A_longi4(1,2)/A_longi4(2,2))
    
    sim('longi.mdl')
    sim('longi3.mdl') 
    sim('longi4.mdl')
    
    figure()
    plot(Temps,incidence,incidence_Ga,incidence_2Gain)
    title('Courbe Incidence \alpha en Vol')
    xlabel('Temps (s)')
    ylabel('Incidence (°)')
    
    figure()
    plot(Temps,incidence,vitesse_tangage_Ga,vitesse_tangage_2Gain)
    title('Courbe Incidence \alpha en Vol')
    xlabel('Temps (s)')
    ylabel('Incidence (°)')
    
    
    