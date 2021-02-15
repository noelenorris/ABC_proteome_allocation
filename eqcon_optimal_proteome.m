%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                          %
%  EQCON_OPTIMAL_PROTEOME.M: Constraints for proteome allocation problem   %
%                                                                          %
%       Specificiation of protoeme allocation model.                       %
%                                                                          %
%       This code is called by function optimal_proteome.m and specifies   %
%       all the constraints for fmincon.                                   %
%                                                                          %
%  Noele Norris                                                            %
%                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [c,ceq]=eqcon_optimal_proteome(x,p,parameter_mat_name)

    % load all specified parameter values
    load(parameter_mat_name)
    
    %external nutrient concentration
    S = x_axis_values(p);
    
    %VARIABLES
    
        %GROWTH RATE
            mu = z(1,p)*x(1);
            
        %CELL RADIUS
            r = z(2,p)*x(2);
            
        %METABOLITES and corresponding molecular weights
        
            %Periplasmic concentration of nutrient
            S_p = z(3,p)*x(3);
            m_S_p = m_Sc;
            
            %Cytoplasmic concentration of nutrient
            S_c = z(4,p)*x(4);

            %Cytoplasmic concentration of amino acids
            A = z(5,p)*x(5);

            %Cytoplasmic concentration of cell wall proteins
            W = z(6,p)*x(6);

            %Cytoplasmic concentration of proteins
            P = z(7,p)*x(7);
            m_P = m_A;
                             
            
            %PROTEOME ALLOCATIONS
            phi_E = z(8,p)*x(8);
            phi_R = z(9,p)*x(9);
            phi_M = z(10,p)*x(10);
            phi_T = z(11,p)*x(11);
            phi_BP = z(12,p)*x(12);
            
            phi = [phi_E; phi_R; phi_M; phi_T; phi_BP];
               
    
        %VOLUME AND SURFACE AREA RATIOS
            vol_ratio_periplasm = z(13,p)*x(13);
            
       
        if(transport_model > 0) %if cell has a periplasm
            
            % ratio of cytoplasmic to periplasmic volume
            ratio_vol_cyto_peri = (1-vol_ratio_periplasm)/vol_ratio_periplasm;
            
            % ratio of volume to surface area of cytoplasm
            ratio_vol_sa_cyto = (1-vol_ratio_periplasm)^(1/n_shape)*r/n_shape;
            
            % ratio of cytoplasmic to periplasmic surface area
            ratio_vol_sa_cyto_peri = (1-vol_ratio_periplasm)*r/n_shape;
        else
            %ratio of volume to surface area of cytoplasm
            ratio_vol_sa_cyto = r/n_shape;
        end
        
               
       %% PROTEOME FRACTIONS 
       
        %CONCENTRATIONS IN THE CYTOPLASM
        
            % catabolic enzymes
            E = P*phi_E/num_aa_E;

            % ribosomes
            R = P*phi_R/num_aa_R;
            
            % cell membrane biosynthesis proteins
            M = P*phi_M/num_aa_M;
            
            % membrane-bound transporters
            T = P*phi_T/num_aa_T;
            
            % binding proteins 
            BP = P*phi_BP/num_aa_BP;
        
        % Convert to periplasmic concentrations
        if(transport_model > 0)
            
            %periplasmic concentration of binding protein
            BP_peri = ratio_vol_cyto_peri*BP;
            
            %periplasmic concenration of membrane-bound transporters
            T_peri = ratio_vol_cyto_peri*T;
        else
            BP_peri = 0;
        end
           
        
    %% Michaelis-Menten enzyme kinetics
        
        % stoichiometry matrix
        N_int = N(1:5,:);
        
        
        %RATES ([cytoplasmic concentration]/[time])
        
            % rate of catabolism
            v1 = k_cat_E*E*S_c/(K_M_E+S_c);
            
            % rate of cell membrane synthesis
            v2 = k_cat_M*M*A/(K_M_M+A);
            
            % rate of protein elongation
            v3 = k_cat_R*R*A/(K_M_R+A);
            
            
            % CYTOPLASMIC TRANSPORT RATE
            
            if(transport_model == 0) %No Periplasm
                
                % Michaelis-Menten constants
                k_cat_T = k2;
                K_M_T = k2/k1;
                
                % cytoplasmic transport rate ([S]_ext --> S_c)
                b1 = k_cat_T*T*S/(K_M_T+S);
                
                S_BP = 0;
                
            elseif(transport_model == 1) % ABC transport system

                BP_total = BP_peri;
                T_total = T_peri;
                k1f = k1;

                % analytical solution of ABC transport model
                T_S_BP = (k3*(k2*k3*k0r - (BP_total^2*S_p^2*k2^2*k0f^2*k1f^2 + 2*BP_total^2*S_p^2*k2*k3*k0f^2*k1f^2 + BP_total^2*S_p^2*k3^2*k0f^2*k1f^2 - 2*BP_total*S_p^2*T_total*k2^2*k0f^2*k1f^2 - 4*BP_total*S_p^2*T_total*k2*k3*k0f^2*k1f^2 - 2*BP_total*S_p^2*T_total*k3^2*k0f^2*k1f^2 + 2*BP_total*S_p^2*k2^2*k3*k0f^2*k1f + 2*BP_total*S_p^2*k2*k3^2*k0f^2*k1f - 2*BP_total*S_p*T_total*k2^2*k3*k0f*k1f^2 - 2*BP_total*S_p*T_total*k2*k3^2*k0f*k1f^2 + 2*BP_total*S_p*k2^2*k3*k0f*k1f*k0r + 2*BP_total*S_p*k2*k3^2*k0f*k1f*k0r + S_p^2*T_total^2*k2^2*k0f^2*k1f^2 + 2*S_p^2*T_total^2*k2*k3*k0f^2*k1f^2 + S_p^2*T_total^2*k3^2*k0f^2*k1f^2 + 2*S_p^2*T_total*k2^2*k3*k0f^2*k1f + 2*S_p^2*T_total*k2*k3^2*k0f^2*k1f + S_p^2*k2^2*k3^2*k0f^2 + 2*S_p*T_total^2*k2^2*k3*k0f*k1f^2 + 2*S_p*T_total^2*k2*k3^2*k0f*k1f^2 + 2*S_p*T_total*k2^2*k3^2*k0f*k1f + 2*S_p*T_total*k2^2*k3*k0f*k1f*k0r + 2*S_p*T_total*k2*k3^2*k0f*k1f*k0r + 2*S_p*k2^2*k3^2*k0f*k0r + T_total^2*k2^2*k3^2*k1f^2 + 2*T_total*k2^2*k3^2*k1f*k0r + k2^2*k3^2*k0r^2)^(1/2) + S_p*k2*k3*k0f + T_total*k2*k3*k1f + BP_total*S_p*k2*k0f*k1f + BP_total*S_p*k3*k0f*k1f + S_p*T_total*k2*k0f*k1f + S_p*T_total*k3*k0f*k1f))/(2*k1f*(k2 + k3)*(k2*k3 + S_p*k2*k0f + S_p*k3*k0f));
                S_BP = -(k2*k3*k0r - (BP_total^2*S_p^2*k2^2*k0f^2*k1f^2 + 2*BP_total^2*S_p^2*k2*k3*k0f^2*k1f^2 + BP_total^2*S_p^2*k3^2*k0f^2*k1f^2 - 2*BP_total*S_p^2*T_total*k2^2*k0f^2*k1f^2 - 4*BP_total*S_p^2*T_total*k2*k3*k0f^2*k1f^2 - 2*BP_total*S_p^2*T_total*k3^2*k0f^2*k1f^2 + 2*BP_total*S_p^2*k2^2*k3*k0f^2*k1f + 2*BP_total*S_p^2*k2*k3^2*k0f^2*k1f - 2*BP_total*S_p*T_total*k2^2*k3*k0f*k1f^2 - 2*BP_total*S_p*T_total*k2*k3^2*k0f*k1f^2 + 2*BP_total*S_p*k2^2*k3*k0f*k1f*k0r + 2*BP_total*S_p*k2*k3^2*k0f*k1f*k0r + S_p^2*T_total^2*k2^2*k0f^2*k1f^2 + 2*S_p^2*T_total^2*k2*k3*k0f^2*k1f^2 + S_p^2*T_total^2*k3^2*k0f^2*k1f^2 + 2*S_p^2*T_total*k2^2*k3*k0f^2*k1f + 2*S_p^2*T_total*k2*k3^2*k0f^2*k1f + S_p^2*k2^2*k3^2*k0f^2 + 2*S_p*T_total^2*k2^2*k3*k0f*k1f^2 + 2*S_p*T_total^2*k2*k3^2*k0f*k1f^2 + 2*S_p*T_total*k2^2*k3^2*k0f*k1f + 2*S_p*T_total*k2^2*k3*k0f*k1f*k0r + 2*S_p*T_total*k2*k3^2*k0f*k1f*k0r + 2*S_p*k2^2*k3^2*k0f*k0r + T_total^2*k2^2*k3^2*k1f^2 + 2*T_total*k2^2*k3^2*k1f*k0r + k2^2*k3^2*k0r^2)^(1/2) + S_p*k2*k3*k0f + T_total*k2*k3*k1f - BP_total*S_p*k2*k0f*k1f - BP_total*S_p*k3*k0f*k1f + S_p*T_total*k2*k0f*k1f + S_p*T_total*k3*k0f*k1f)/(2*k1f*(k2 + k3)*(k0r + S_p*k0f));

                 % cytoplasmic transport rate ([S]_p --> S_c)
                 b1 = k2*T_S_BP/ratio_vol_cyto_peri;
                
                
            else %PTS 
                
                %Micahelis-Menten constants
                k_cat_T = k2;
                K_M_T = k2/k1;
                
                % cytoplasmic transport rate ([S]_p --> S_c)
                b1 = k_cat_T*T*S_p/(K_M_T+S_p);
                
                S_BP = 0;
            end
            
            
            % PERIPLASMIC DIFFUSIVE UPTAKE RATE ([S] --> [S]_p)
            if(transport_model == 0)
                
                b2 = 0;
            else
                
                % periplasmic uptake rate [cytoplasmic concentration/time]
                % NB: assuming sphere to calculate diffusive flux rate
                b2 = 3*D_substrate*(S-S_p)*r^(-2)/(1-vol_ratio_periplasm);
                
            end

            % vector of enzymatic rates
            v_rates = [v1; v2; v3; b1; b2];
        
            
        %% METABOLITES
        
        concentrations = [S_c; A; W; P; S_p/ratio_vol_cyto_peri];
        
        % fraction of proteome freely diffusing in cytoplasm
        phi_cyto = phi_O_cyto+phi_E+phi_R+phi_M;
        
        % fraction of proteome freely diffusing in periplasm
        phi_peri = phi_O_peri+phi_BP;
      
        % cytoplamic concentrations
        concentrations_cyto = [S_c; A; W; P*phi_cyto];
        mol_weights_cyto = [m_Sc; m_A; m_w; m_P];
        
        % periplasmic concentrations
        if(transport_model > 0)
            concentrations_peri = [P*phi_peri*ratio_vol_cyto_peri; S_p+S_BP];
            mol_weights_peri = [m_P; m_S_p];
        end


        %% STEADY-STATE BALANCED,GROWTH-RATE CONSTRAINTS
        ceq(1:5) = 1000*(N_int*v_rates - mu*concentrations);
        
    
        %% PROTEOME ALLOCATION CONSTRAINT
        ceq(6) = 1 - phi_O_cyto-phi_O_peri-phi_O_IM - sum(phi);

    
        %% CELL WALL CONSTRAINT
        if(transport_model == 0) %if no periplasm, only need inner membrane
            ceq(7) = 1 - W*a_w/(1/ratio_vol_sa_cyto);
        else %if periplasm, need cell membrane units for both inner and outer membranes
            ceq(7) = 1 - W*a_w/(1/ratio_vol_sa_cyto+1/ratio_vol_sa_cyto_peri);
        end
    
    
        %% CYTOPLASMIC DENSITY CONSTRAINT
        c(1) = concentrations_cyto'*mol_weights_cyto - density_available_cyto*density;
    
        %% MEMBRANE-BOUND TRANSPORT UNIT REAL_ESTATE CONSTRAINT
        c(2) = ratio_vol_sa_cyto*(T+P*phi_O_IM/num_aa_T)*a_transporter-SA_available_cyto;
    

        if(transport_model > 0)
            %% PERIPLASMIC DENSITY CONSTRAINT
            c(3) = concentrations_peri'*mol_weights_peri - density_available_peri*density;
        end
end
 