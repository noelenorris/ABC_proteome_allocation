clear all;
close all;
clc;

seed = rng('shuffle');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                            %
%  MAIN.M: A PROTEOME ALLOCATION PROBLEM WITH ABC OR PTS TRANSPORT SYSTEMS   %
%                                                                            %
%                                                                            %
%       !IMPORTANT! PLEASE READ: This code will not run  without             %
%       MATLAB's Global Optimization Toolbox. To install the toolbox,        %
%       click: Home --> Add-Ons --> Get Add-Ons.                             %
%       Search for "Global Optimization Toolbox" and click on toolbox        %
%       then "Install". Follow the prompts to finish installation.           %
%       (In older versions of MATLAB, you can obtain the toolbox             %
%        by re-running the installer.)                                       %
%                                                                            %
%       This file specifies all parameter values used to solve a proteome    %
%       allocation problem. It saves the parameter values in a .mat file     %
%       and calls the optimal_proteome function to run the allocation        %
%       problem.                                                             %
%                                                                            %
%  Noele Norris                                                              %
%                                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% NB: PARAMETER UNITS
%
%     Concentration: umoles / mL  = 1 mM (mL = cm^3)
%     Length scale: cm
%     Time scale: msec
%     Mass: decigram
%
%%
    
 
for j = [0.1 10]   %(for sensitivity analysis, iterate over parameter values)
    
    %% OPTIONS 
        % cell shape for calculating volume and surface area ratios
        n_shape = 3; %3: sphere, 2: infinitely long cylinder

        % specify transport system of cell
        transport_model = 1;  %0: No Periplasm, 1: ABC, 2: PTS

        %file name of .mat which stores parameters specified in this file
        mat_file_name = strrep(strcat('ABC_BP_cost_', num2str(j)),'.','p'); 
        
        % number of iterations to repeat optimization with different IC
        num_it = 50; 
    

    
    %% CELL PARAMETER VALUES 
    % See supplemental information for justification of chosen parameter
    % values.
   
    
        %AVAILABLE CELL DENSITY

            %basal proteomic requirement for different cell compartments
              phi_O_cyto = 0;
              phi_O_peri = 0;
              phi_O_IM = 0;

            %dry mass cell density
            density = 0.34*10^1; % decigrams/mL (BioNumbers ID: 109049)
            density_available_cyto = 1; % fraction of density available for cyptolasm
            density_available_peri = 1; % fraction of density available for periplasm



        %AVAILABLE SURFACE AREA
            SA_available_cyto = 0.5; %fraction of inner membrane surface available for proteins


        % ESTIMATED PROTEOMIC COSTS

            num_aa_T = 1600; %membrane-bound transport unit
            num_aa_BP = j*400; %binding protein
            
            num_aa_E = 40000; % metabolism protein unit
            num_aa_M = 21000; % cell wall synthesis protein unit
            num_aa_R = 20000; %protein synthesis unit



        % MOLECULAR WEIGHTS (dag per umole)
        
            m_A = 110*10^-6*10^1; %amino acid
            m_Sc = 180.1*10^-6*10^1; %intracellular metabolite 
            m_w = 1*220*10^1*10^-6; %cell membrane unit

            Av = 6.0221409e+23*10^-6; % Avogadro's number (numbers/umole)


        % REQUIRED SURFACE AREA (cm^2/umol)

            a_transporter = 20*10^-14*Av; %membrane-bound transport unit
            a_w = 2.5e-14*Av; %cell membrane unit


        % MICHAELIS-MENTEN KINETICS PARAMETER VALUES


            % BINDING PROTEIN
            %K_BP = (0.1*0.001);
            k0f = 10^5/1000;
            k0r = 100/1000;

            % CYTOPLASMIC TRANSPORT
            % V_max = 100 umoles/min/mg
            if(transport_model == 1) % ABC transporter
                 k2 = (210/1000); %msec^-1
                 k1 = .01*21;
                 k3 = (210/1000)*0.01;
            else % PTS transporter
                 k2 = (210/1000); 
                 k1 = 21;
            end

            % CELL WALL BIOSYNTHESIS  
            k_cat_M = 3.14/1000; %msec^-1
            K_M_M = 20*0.001;  %umoles/mL

            % METABOLISM (E): S_c --> Amino acid
            k_cat_E = 100/1000; % msec^-1
            K_M_E = 58*0.001; %umol/mL 

            % PROTEIN SYNTHESIS
            k_cat_R = 15/1000; % msec^-1
            K_M_R = 300*.001; %umol/mL



        % REACTION STOICHIOMETRIES 

            %     v1  v2  v3    b1   b2   
             N = [-5   0   0   1    0 ; ... % S_c   
                   6  -1  -1   0    0 ; ... %  A
                   0   1   0   0    0 ; ... %  W
                   0   0   1   0    0 ; ... %  P
                   0   0   0  -1    1 ;  ... % S_p
                   0   0   0   0   -1];      %S

            %S --> S_c (if there is no periplasm)
            if(transport_model == 0)
                N(5,4) = 0;
                N(6,4) = -1;
            end
                       
        
    %% ENVIRONMENT PARAMETER VALUES
  
        x_axis_values = 10.^[-6:1:2]; %range of 1nM to 100mM 
        
        D_substrate = 600*10^(-8); %cm^2/sec

        
    %% LOWER AND UPPER BOUNDS for optimization constraints
        
            if(transport_model==1)
                r_max = 0.0001;
            else
               r_max = 0.0001*100;
            end
            
            % Constrain radius to a minimum of 60 nm.
            r_min = 0.0001*0.06;

            % At least one metabolite within cell
            r_met = 0.0001;
            lb_metabolite = 1/(4*pi*r_met^3*Av/3);
            ub_metabolite = 10^6;
            
            
            if(transport_model == 0)
                min_bounds = [0; r_min; 0; lb_metabolite*ones(4,1); 0*ones(5,1); 0];
            else
                min_bounds = [0; r_min; lb_metabolite*ones(5,1); 0*ones(5,1); 0];
            end

           
            max_bounds = [1000; r_max; ub_metabolite*ones(5,1); (1-phi_O_cyto-phi_O_IM-phi_O_peri)*ones(5,1); 1];

            
            if(transport_model == 0) %if there is no periplasm
                max_bounds(3) = 0;
                max_bounds(12) = 0;
                min_bounds(13) = 0;
                max_bounds(13) = 0;
            end
            
            if(transport_model == 2)
                max_bounds(12) = 0;
            end
            
            
        
   %% VARIABLE TRANSFORMATION
   % Transform variables based on predicted optimal solutions for
   % fmincon tractability.

        z_1 = ones(13,1);   

        z_1(1) = 10^(-6);   %(1) mu
        z_1(2) = 10^(-4);   %(2)  r
        z_1(4) = 10^1;      %(4) [S_c]    
        z_1(5) = 10^1;      %(5) [A] 
        z_1(6) = 10^1;      %(6) [W]
        z_1(7) = 10^3;      %(7) [P]
        z_1(8) = 10^(-2);   %(8) metabolism, phi_E
        z_1(9) = 10^(0);    %(9) protein synthesis, phi_R
        z_1(10) = 10^(-1);  %(10) cell wall synthesis, phi_M
        z_1(11) = 10^(0);   %(11) membrane-bound transport, phi_T
        z_1(12) = 10^(-1);  %(12) binding protein, phi_BP
        z_1(13) = 10^(0);   %(13) periplasmic fraction of volume, f_peri

        z = repmat(z_1, [1 length(x_axis_values)]);

        z(3,:) = x_axis_values; %(3) [S_p]
 
     
  
   %%
      save(mat_file_name); % save parameter values in .mat
      
      optimal_proteome(mat_file_name); % run optimization using saved parameters
      
end
          
