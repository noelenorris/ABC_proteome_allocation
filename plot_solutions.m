
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                          %
%  PLOT_SOLUTIONS.M: Create various figures comparing solutions            %
%                                                                          %
% of proteome allocation problem with ABC transport or PTS.                %
%                                                                          %
%                                                                          %
%       Noele Norris                                                       %
%                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


main_title_k{1} = 'rho_peri'
fig_filename_k{1} = 'rho_peri'
results_ABC{1} = 'Results/results_January152021120903701_ABC_phi_peri_x0p01'
results_ABC_2{1} = 'Results/results_January152021134746540_ABC_phi_peri_x100'
labels_k{1} = {'PTS', 'Baseline', '\rho_{peri} x0.01', '\rho_{peri} x100'};

main_title_k{2} = 'rho_cyto'
fig_filename_k{2} = 'rho_cyto'
results_ABC{2} = 'Results/results_January152021092050607_ABC_phi_cyto_x0p01'
results_ABC_2{2} = 'Results/results_January152021101446536_ABC_phi_cyto_x100'
labels_k{2} = {'PTS', 'Baseline', '\rho_{cyto} x0.01', '\rho_{cyto} x100'};




for k=1:1:length(labels_k)

% % BASELINE PLOTS
% 
%     % specify which results to plot
%         results = { ...
%         'Results/results_January122021081212150_PTS1', ...
%         'Results/results_January112021162212132_ABC', ....
%     };
%     labels = {'PTS', 'ABC'};
%     main_title = 'Baseline';
%     fig_filename = 'Baseline';

  % SENSITIVITY ANALYSIS
    results = { ...
        'Results/results_January122021081212150_PTS1', ...
        'Results/results_January112021162212132_ABC', ...
        results_ABC{k}, ...
        results_ABC_2{k}, ...
    };

    labels = labels_k{k};
    main_title = main_title_k{k};
    fig_filename = fig_filename_k{k};



    % colormap
    col = lines(7);


   %% Process results       

    num_results = length(results);
    for it=1:1:num_results

        load(results{it})

        opt_values = cell2mat(optimal_values);

        opt_val_array{it} = opt_values;


        z_val = z';

        S{it} = x_axis_values;

        mu{it} = z_val(:,1).*opt_values(:,1);
        r{it} = z_val(:,2).*opt_values(:,2);
        S_p{it} = z_val(:,3).*opt_values(:,3);
        S_c{it} = z_val(:,4).*opt_values(:,4);
        A{it} = z_val(:,5).*opt_values(:,5);
        W{it} = z_val(:,6).*opt_values(:,6);
        P{it} = z_val(:,7).*opt_values(:,7);
        phi_E{it} = z_val(:,8).*opt_values(:,8);
        phi_R{it} = z_val(:,9).*opt_values(:,9);
        phi_M{it} = z_val(:,10).*opt_values(:,10);
        phi_T{it} = z_val(:,11).*opt_values(:,11);
        phi_BP{it} = z_val(:,12).*opt_values(:,12);
        vol_ratio_periplasm{it} = z_val(:,13).*opt_values(:,13);

        proteome_fraction_sums = phi_E{it}+phi_R{it}+phi_M{it}+phi_T{it}+phi_BP{it}

        % Calculate Michaelis-Menten transport parameters
            ratio_vol_cyto_peri = (1-vol_ratio_periplasm{it})./vol_ratio_periplasm{it};


           %CONCENTRATIONS IN THE CYTOPLASM
            E = P{it}.*phi_E{it}/num_aa_E;
            R = P{it}.*phi_R{it}/num_aa_R;
            M = P{it}.*phi_M{it}/num_aa_M;
            T = P{it}.*phi_T{it}/num_aa_T;
            BP = P{it}.*phi_BP{it}/num_aa_BP;

            ratio_num_BP_num_T{it} = BP./T;

            %CONCENTRATIONS IN THE PERIPLASM
            if(transport_model > 0)
                BP_peri = ratio_vol_cyto_peri.*BP;
                T_peri = ratio_vol_cyto_peri.*T;
            else
                BP_peri = 0;
            end

            % STOICHIOMETRIC MATRIX   
            N_int = N(1:5,:);

            % RATES
                v1 = k_cat_E.*E.*S_c{it}./(K_M_E+S_c{it});
                v2 = k_cat_M.*M.*A{it}./(K_M_M+A{it});
                v3 = k_cat_R.*R.*A{it}./(K_M_R+A{it});
                
            
           % CALCULATE UPTAKE RATE
            %Calculate V_max
              S_p2 = S_p{it};
              if(transport_model == 1) %ABC transport model

                    BP_total = BP_peri;
                    T_total = T_peri;
                    k1f = k1;

                    T_S_BP = (k3.*(k2.*k3.*k0r - (BP_total.^2.*S_p2.^2.*k2.^2.*k0f.^2.*k1f.^2 + 2.*BP_total.^2.*S_p2.^2.*k2.*k3.*k0f.^2.*k1f.^2 + BP_total.^2.*S_p2.^2.*k3^2*k0f^2*k1f^2 ...
                           - 2*BP_total.*S_p2.^2.*T_total.*k2^2.*k0f^2.*k1f^2 - 4.*BP_total.*S_p2.^2.*T_total*k2*k3*k0f^2*k1f^2 - 2*BP_total.*S_p2.^2.*T_total.*k3^2*k0f^2*k1f^2 ...
                           + 2*BP_total.*S_p2.^2.*k2^2*k3*k0f^2*k1f + 2*BP_total.*S_p2.^2.*k2*k3^2*k0f^2*k1f - 2*BP_total.*S_p2.*T_total.*k2^2*k3*k0f*k1f^2 - 2*BP_total.*S_p2.*T_total.*k2*k3^2*k0f*k1f^2 ...
                           + 2*BP_total.*S_p2.*k2^2*k3*k0f*k1f*k0r + 2*BP_total.*S_p2*k2*k3^2*k0f*k1f*k0r + S_p2.^2.*T_total.^2*k2^2*k0f^2*k1f^2 + 2*S_p2.^2.*T_total.^2*k2*k3*k0f^2*k1f^2 ...
                           + S_p2.^2.*T_total.^2*k3^2*k0f^2*k1f^2 + 2.*S_p2.^2.*T_total*k2^2*k3*k0f^2*k1f + 2*S_p2.^2.*T_total*k2*k3^2*k0f^2*k1f + S_p2.^2*k2^2*k3^2*k0f^2 + 2*S_p2.*T_total.^2*k2^2*k3*k0f*k1f^2 ...
                           + 2*S_p2.*T_total.^2*k2*k3^2*k0f*k1f^2 + 2*S_p2.*T_total*k2^2*k3^2*k0f*k1f + 2*S_p2.*T_total*k2^2*k3*k0f*k1f*k0r + 2*S_p2.*T_total*k2*k3^2*k0f*k1f*k0r + 2*S_p2.*k2^2*k3^2*k0f*k0r ...
                           + T_total.^2*k2^2*k3^2*k1f^2 + 2*T_total.*k2^2*k3^2*k1f*k0r + k2^2*k3^2*k0r^2).^(1/2) + S_p2.*k2*k3*k0f + T_total.*k2*k3*k1f + BP_total.*S_p2*k2*k0f*k1f + BP_total.*S_p2*k3*k0f*k1f ...
                           + S_p2.*T_total*k2*k0f*k1f + S_p2.*T_total*k3*k0f*k1f))./(2*k1f*(k2 + k3)*(k2*k3 + S_p2.*k2*k0f + S_p2.*k3*k0f));

                     % cytoplasmic transport rate ([S]_p --> S_c)
                     uptake{it} = k2*T_S_BP./ratio_vol_cyto_peri;

            else % PTS
                        
                    k_cat_T = k2;
                    
                    %[S]_p --> S_c
                    uptake{it} = k_cat_T*T.*S_p2./(k2/k1+S_p2)
                    S_BP = 0;
            end


              % CALCULATE V_MAX
              S_p2 = 10^10;
              if(transport_model == 1) %ABC 

                    %% ABC Transport Model

                    BP_total = BP_peri;
                    T_total = T_peri;
                    k1f = k1;

                    % NEW SOLUTION
                    T_S_BP = (k3.*(k2.*k3.*k0r - (BP_total.^2.*S_p2.^2.*k2.^2.*k0f.^2.*k1f.^2 + 2.*BP_total.^2.*S_p2.^2.*k2.*k3.*k0f.^2.*k1f.^2 + BP_total.^2.*S_p2.^2.*k3^2*k0f^2*k1f^2 ...
                           - 2*BP_total.*S_p2.^2.*T_total.*k2^2.*k0f^2.*k1f^2 - 4.*BP_total.*S_p2^2.*T_total*k2*k3*k0f^2*k1f^2 - 2*BP_total.*S_p2.^2.*T_total.*k3^2*k0f^2*k1f^2 ...
                           + 2*BP_total.*S_p2.^2.*k2^2*k3*k0f^2*k1f + 2*BP_total.*S_p2.^2.*k2*k3^2*k0f^2*k1f - 2*BP_total.*S_p2.*T_total.*k2^2*k3*k0f*k1f^2 - 2*BP_total.*S_p2.*T_total.*k2*k3^2*k0f*k1f^2 ...
                           + 2*BP_total.*S_p2.*k2^2*k3*k0f*k1f*k0r + 2*BP_total.*S_p2*k2*k3^2*k0f*k1f*k0r + S_p2.^2.*T_total.^2*k2^2*k0f^2*k1f^2 + 2*S_p2.^2.*T_total.^2*k2*k3*k0f^2*k1f^2 ...
                           + S_p2.^2.*T_total.^2*k3^2*k0f^2*k1f^2 + 2.*S_p2.^2.*T_total*k2^2*k3*k0f^2*k1f + 2*S_p2^2*T_total*k2*k3^2*k0f^2*k1f + S_p2^2*k2^2*k3^2*k0f^2 + 2*S_p2.*T_total.^2*k2^2*k3*k0f*k1f^2 ...
                           + 2*S_p2.*T_total.^2*k2*k3^2*k0f*k1f^2 + 2*S_p2.*T_total*k2^2*k3^2*k0f*k1f + 2*S_p2.*T_total*k2^2*k3*k0f*k1f*k0r + 2*S_p2.*T_total*k2*k3^2*k0f*k1f*k0r + 2*S_p2.*k2^2*k3^2*k0f*k0r ...
                           + T_total.^2*k2^2*k3^2*k1f^2 + 2*T_total.*k2^2*k3^2*k1f*k0r + k2^2*k3^2*k0r^2).^(1/2) + S_p2.*k2*k3*k0f + T_total.*k2*k3*k1f + BP_total.*S_p2*k2*k0f*k1f + BP_total.*S_p2*k3*k0f*k1f ...
                           + S_p2.*T_total*k2*k0f*k1f + S_p2.*T_total*k3*k0f*k1f))./(2*k1f*(k2 + k3)*(k2*k3 + S_p2.*k2*k0f + S_p2.*k3*k0f));


                     % cytoplasmic transport rate ([S]_p --> S_c)
                     b1 = k2*T_S_BP./ratio_vol_cyto_peri;

                else %PTS
                    
                    k_cat_T = k2;
                    %[S]_p --> S_c
                    b1 = k_cat_T*T.*S_p2./(k2/k1+S_p2);
                    S_BP = 0;
              end


               V_max_T{it} = b1; %mM/msec

               V_max_downstream{it} = -N(1,1)*k_cat_E.*E;

               % Calculate K_effective
               if(transport_model == 1)

                   S_p2 = 10.^[-10:.01:-2];
                   for it2=1:1:length(BP_peri)
                        BP_total = BP_peri(it2);
                        T_total = T_peri(it2);
                        k1f = k1;

                        T_S_BP = (k3.*(k2.*k3.*k0r - (BP_total.^2.*S_p2.^2.*k2.^2.*k0f.^2.*k1f.^2 + 2.*BP_total.^2.*S_p2.^2.*k2.*k3.*k0f.^2.*k1f.^2 + BP_total.^2.*S_p2.^2.*k3^2*k0f^2*k1f^2 ...
                       - 2*BP_total.*S_p2.^2.*T_total.*k2^2.*k0f^2.*k1f^2 - 4.*BP_total.*S_p2.^2.*T_total*k2*k3*k0f^2*k1f^2 - 2*BP_total.*S_p2.^2.*T_total.*k3^2*k0f^2*k1f^2 ...
                       + 2*BP_total.*S_p2.^2.*k2^2*k3*k0f^2*k1f + 2*BP_total.*S_p2.^2.*k2*k3^2*k0f^2*k1f - 2*BP_total.*S_p2.*T_total.*k2^2*k3*k0f*k1f^2 - 2*BP_total.*S_p2.*T_total.*k2*k3^2*k0f*k1f^2 ...
                       + 2*BP_total.*S_p2.*k2^2*k3*k0f*k1f*k0r + 2*BP_total.*S_p2*k2*k3^2*k0f*k1f*k0r + S_p2.^2.*T_total.^2*k2^2*k0f^2*k1f^2 + 2*S_p2.^2.*T_total.^2*k2*k3*k0f^2*k1f^2 ...
                       + S_p2.^2.*T_total.^2*k3^2*k0f^2*k1f^2 + 2.*S_p2.^2.*T_total*k2^2*k3*k0f^2*k1f + 2*S_p2.^2.*T_total*k2*k3^2*k0f^2*k1f + S_p2.^2*k2^2*k3^2*k0f^2 + 2*S_p2.*T_total.^2*k2^2*k3*k0f*k1f^2 ...
                       + 2*S_p2.*T_total.^2*k2*k3^2*k0f*k1f^2 + 2*S_p2.*T_total*k2^2*k3^2*k0f*k1f + 2*S_p2.*T_total*k2^2*k3*k0f*k1f*k0r + 2*S_p2.*T_total*k2*k3^2*k0f*k1f*k0r + 2*S_p2.*k2^2*k3^2*k0f*k0r ...
                       + T_total.^2*k2^2*k3^2*k1f^2 + 2*T_total.*k2^2*k3^2*k1f*k0r + k2^2*k3^2*k0r^2).^(1/2) + S_p2.*k2*k3*k0f + T_total.*k2*k3*k1f + BP_total.*S_p2*k2*k0f*k1f + BP_total.*S_p2*k3*k0f*k1f ...

                        uptake_K_eff = k2*T_S_BP./ratio_vol_cyto_peri(it2);

                        if(~isnan(uptake_K_eff(1)))

                            f_uptake = fit(S_p2', (uptake_K_eff/V_max_T{it}(it2))', 'x/(a+x)', 'Lower', [0], 'Upper', [10^(-3)], 'StartPoint', [10^(-5)]);

                            K_eff_it(it2) = f_uptake.a;
                        else
                            K_eff_it(it2) = NaN;
                        end

                   end

                   K_eff{it} = K_eff_it';


               else
                   K_eff{it} = k2/k1*ones(size(b1));
               end



    end


    %%

% % 
%     % %% Maximal growth rate vs. external nutrient concentration
%         figure()
%         it = 1;
%         semilogx(S{it}, 3.6*10^6*mu{it},'-o','LineWidth',num_results+5-it) 
%         hold on
%         while(it < num_results)
%              it = it + 1;
%              loglog(S{it}, 3.6*10^6*mu{it},'-o','LineWidth', num_results + 5-it)
%         end
%         xlabel('Nutrient concentration [mM]')
%         xlim([10^(-6) 10^2])
%         xticks([10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
%         ylabel('Maximal growth rate [1/hour]')
%         set(findall(gcf,'-property','FontSize'),'FontSize',25)
%         grid on 
%         hold off
%     
%     
%     % %% Optimal radius vs. external nutrient concentration
%     %     figure(2)
%     %     it = 1;
%     %     loglog(S{it}, 10^4*r{it},'-o','LineWidth',num_results+5-it)
%     %     hold on
%     %     while(it < num_results)
%     %         it = it + 1;
%     %         loglog(S{it}, 10^4*r{it},'-o','LineWidth',num_results+5-it)
%     %     end
%     %     xlabel('Nutrient concentration [mM]')
%     %     xlim([10^(-6) 10^2])
%     %     xticks([10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
%     %     ylabel('Optimal cell radius (\mum)')
%     %     legend(labels, 'Location', 'southeast')   
%     %     set(findall(gcf,'-property','FontSize'),'FontSize',25)
%     %     grid on 
%     %     hold off
%     
%         
%         %% Optimal SA/Vol vs. external nutrient concentration
%         figure()
%         it = 1;
%         loglog(S{it}, 3./(10^4*r{it}),'-o','LineWidth',num_results+5-it)
%         hold on
%         while(it < num_results)
%             it = it + 1;
%             loglog(S{it}, 3./(10^4*r{it}),'-o','LineWidth',num_results+5-it)
%         end
%         xlabel('Nutrient concentration [mM]')
%         xlim([10^(-6) 10^2])
%         xticks([10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
%         ylabel('Optimal SA-to-Volume ratio [1/\mum]')
%         legend(labels, 'Location', 'southeast')   
%         set(findall(gcf,'-property','FontSize'),'FontSize',25)
%         grid on 
%         hold off
%         
%     %% Optimal periplasmic volume fraction vs. external nutrient concentration
%         figure()
%         it = 1;
%         semilogx(S{it}, vol_ratio_periplasm{it},'-o','LineWidth',num_results+5-it)
%         hold on
%         while(it < num_results)
%             it = it + 1;
%             semilogx(S{it}, vol_ratio_periplasm{it},'-o','LineWidth',num_results+5-it)
%         end
%         xlabel('Nutrient concentration [mM]')
%         ylabel('Optimal periplasmic volume fraction', 'FontSize', 16)
%         xticks([10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
%         legend(labels)
%         set(findall(gcf,'-property','FontSize'),'FontSize',25)
%         grid on 
%         hold off
%     
%         
%         
%     %% Optimal proteome allocation vs. external nutrient concentration
%         col=lines(6);
%         fig = figure();
%         ax = axes('parent',fig);
%         loglog(S{1}, phi_T{1}+phi_BP{1},'-','LineWidth',8, 'Color',[0.5843    0.8157    0.9882])
%         hold on
%         loglog(S{2}, phi_T{2}+phi_BP{2},'-','LineWidth',8, 'Color',[.98 .6 .51])
%         loglog(S{1}, phi_T{1},'-','LineWidth',4,'Color',col(1,:))
%         loglog(S{2}, phi_T{2},'-','LineWidth',4, 'Color',col(2,:))
%         loglog(S{2}, phi_BP{2},'--','LineWidth',3, 'Color',col(2,:))
%         loglog(S{1}, phi_E{1},':','LineWidth',4, 'Color',col(1,:))
%         loglog(S{2}, phi_E{2},':','LineWidth',4, 'Color',col(2,:))
%         loglog(S{1}, phi_R{1}, '--','LineWidth',2, 'Color', col(1,:))
%         loglog(S{it}, phi_R{2}, '--','LineWidth',2, 'Color', col(2,:))
%         xlim([10^-6 10^2])
%         xticks([10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
%         xlabel('Nutrient concentration [mM]', 'FontSize', 16)
%         ylabel('Proteome fraction','FontSize', 16)
%         title('Optimal Proteome Allocation', 'FontSize', 16)
%         dummyax = axes('parent',gcf,'position',get(ax,'position'));
%         hold(dummyax, 'on');
%         set(dummyax,'visible','off')
%         h1 = plot(dummyax, nan(1,2), 'LineWidth',3, 'Color', col(1,:));
%         h2 = plot(dummyax, nan(1,2), 'LineWidth',3, 'Color', col(2,:));
%         legend([h1 h2],'PTS','ABC','Location', 'southeast', 'FontSize', 14)
%         dummyax2 = axes('parent',gcf,'position',get(ax,'position'));
%         hold(dummyax2, 'on');
%         set(dummyax2,'visible','off')
%         h7 = plot(dummyax2, nan(1,2), '-', 'LineWidth', 8, 'Color', [0.5 0.5 0.5]);
%         h4 = plot(dummyax2, nan(1,2), '-','LineWidth',4, 'Color', [0 0 0]);
%         h5 = plot(dummyax2, nan(1,2), '--','LineWidth',3, 'Color', [0 0 0]);
%         h6 = plot(dummyax2, nan(1,2), ':','LineWidth',4, 'Color', [0 0 0]);
%         h8 = plot(dummyax2, nan(1,2), '--','LineWidth',2, 'Color', [0 0 0])
%         set(findall(gcf,'-property','FontSize'),'FontSize',25)
%         legend([h7 h4 h5 h6 h8],'Total transport', 'Transport unit', 'Binding protein', 'Metabolism', 'Ribosomes','Location','south', 'FontSize', 14)
%         grid on 
%         hold off
%   
%         
%    %% uptake rate per proteomic cost
%    figure()
%            it = 1;
%         loglog(S{it}, uptake{it}./(phi_T{it}+phi_BP{it}),'-o','LineWidth',num_results+5-it) 
%         hold on
%         while(it < num_results)
%              it = it + 1;
%              loglog(S{it}, uptake{it}./(phi_T{it}+phi_BP{it}),'-o','LineWidth', num_results + 5-it)
%         end
%         xlabel('Nutrient concentration [mM]')
%         xlim([10^(-6) 10^2])
%         xticks([10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
%         ylabel('Uptake rate per transport proteome fraction [mM/hour]')
%         legend(labels)
%         set(findall(gcf,'-property','FontSize'),'FontSize',25)
%         grid on 
%         hold off
%    
%             
%         
%         
%     
%     %% Optimal maximal uptake rate and specific affinity vs. external nutrient concentration
%         figure()
%             subplot(2,1,1)
%             it = 1;
%             loglog(S{it}, 3.6*10^6*V_max_T{it},'-o', 'LineWidth',num_results+5-it)
%             hold on
%             while(it < num_results)
%                 it = it + 1;
%                 loglog(S{it}, 3.6*10^6*V_max_T{it},'-o', 'LineWidth',num_results+5-it)
%             end  
%             xlabel('Nutrient concentration [mM]')
%             xlim([10^(-6) 10^2])
%             xticks([10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
%             yticks([10^3 10^4 10^5 10^6 10^7 10^8])
%             ylabel('V_{max} [mM/hour]')
%             title('Optimal maximal uptake rate')
%             legend(labels, 'FontSize', 14)
%             set(findall(gcf,'-property','FontSize'),'FontSize',25)
%             grid on
%             hold off
%             
%            
%         
%             subplot(2,1,2)
%             it = 1;
%             loglog(S{it},3.6*10^6*V_max_T{it}./K_eff{it},'-o', 'LineWidth',num_results+5-it)
%             hold on
%             while(it < num_results)
%                 it = it + 1;
%                 loglog(S{it}, 3.6*10^6*V_max_T{it}./K_eff{it},'-o', 'LineWidth',num_results+5-it)
%             end  
%             xlabel('Nutrient concentration [mM]', 'FontSize', 16)
%             xlim([10^(-6) 10^2])
%             xticks([10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
%             yticks([10^5 10^6 10^7 10^8 10^9 10^10])
%             ylabel('Specific affinity [1/hour]', 'FontSize', 16)
%             legend(labels, 'FontSize', 14, 'Location', 'east')
%             title('Optimal specific affinity', 'FontSize', 16)
% 
%             set(findall(gcf,'-property','FontSize'),'FontSize',25)
%             grid on
%             hold off
%             
%             
%     
%     
%     %% Optimal ratio of maximal uptake and maximal catabolic rate vs. external nutrient concentration
%         figure()
%         it = 1;
%         loglog(S{it},V_max_T{it}./V_max_downstream{it},'-o', 'LineWidth',num_results+5-it)
%         hold on
%         while(it < num_results)
%             it = it + 1;
%             loglog(S{it}, V_max_T{it}./V_max_downstream{it},'-o', 'LineWidth',num_results+5-it)
%         end  
%         xlabel('Nutrient concentration [mM]')
%         xlim([10^(-6) 10^2])
%         xticks([10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
%         ylim([10^(-1) 10^(5)])
%         yticks([10^(-1) 10^(1) 10^(3) 10^(5)])
%         ylabel('Ratio')
%         title('Ratio of Optimal Maximal Uptake and Maximal Catabolic Rates')
%         legend(labels)
%         
%         set(findall(gcf,'-property','FontSize'),'FontSize',25) 
%         grid on
%         hold off
%     
%     %% Half-saturation constant and BP/T ratio 
%         figure()
%         
%         % Colors: [0.64,0.08,0.18], [0.30,0.75,0.93]
%     
%         yyaxis left
%         grid on
%         it = 2;
%         semilogx(S{it},10^6*K_eff{it},'-o', 'LineWidth',num_results+5-it)
%     
%         xlabel('Nutrient concentration [mM]', 'FontSize', 16)
%         xticks([10^(-6) 10^(-4) 10^(-2) 10^(0) 10^(2)])
%         ylabel('K_{eff} [nM] ', 'FontSize', 16)
%         %yticks([1 10 100 1000 1e4 1e5])
%         ylim([0 25])
%         yticks([0 5 10 15 20 25])
%         %title('Optimal Half-Saturation Constant', 'FontSize', 16)
%         
%        
%         
%         
%         yyaxis right
%         semilogx(S{it}, ratio_num_BP_num_T{it},'-o', 'LineWidth',num_results+5-it)
%         ylabel('Ratio of binding proteins to transport units')
%         ylim([0 10])
%         yticks([0 2 4 6 8 10])
%         xticks([10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
%         
%         ax = gca;
%         ax.ColorOrderIndex = 1;
%         dummyax = axes('parent',gcf,'position',get(ax,'position'));
%         hold(dummyax, 'on');
%         set(dummyax,'visible','off')
%         h2 = plot(dummyax, nan(1,2), 'LineWidth',num_results+2, 'Color', [0 0 0]);
%         h1 = plot(dummyax, nan(1,2), ':','LineWidth',num_results+2, 'Color', [0 0 0]);
%         legend([h2],labels{[2]},'Location', 'northwest')
%         
%         set(findall(gcf,'-property','FontSize'),'FontSize',25)
%         
%         hold off

     %%    
        fig = figure()
        subplot(2,2,1)
            it = 1;
            semilogx(S{it}, 3.6*10^6*mu{it},'-o','LineWidth',6-it, 'Color', col(it,:)) 
            hold on
            while(it < num_results)
                 it = it + 1;
                 loglog(S{it}, 3.6*10^6*mu{it},'-o','LineWidth', 6-it, 'Color', col(it,:))
            end
            xlabel('Nutrient concentration [mM]')
            xlim([10^(-6) 10^2])
            xticks([10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
            ylabel('Growth rate [1/hour]')
            title('Maximal growth rate')
            set(findall(gca,'-property','FontSize'),'FontSize',16)
            legend(labels(1:end), 'Location', 'southeast', 'FontSize', 14)
            grid on 
            hold off
        subplot(2,2,2)
            it = 2;
            semilogx(S{it}, 3./(10^4*r{it}),'-o','LineWidth',6-it, 'Color', col(it,:))
            hold on
            while(it < num_results)
                it = it + 1;
                semilogx(S{it}, 3./(10^4*r{it}),'-o','LineWidth',6-it, 'Color', col(it,:))
            end
            xlabel('Nutrient concentration [mM]')
            xlim([10^(-6) 10^2])
            xticks([10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
            title('Optimal SA-to-volume ratio')
            ylabel('Ratio [1/\mum]')
            set(findall(gca,'-property','FontSize'),'FontSize',16)
            legend(labels(2:end), 'Location', 'northeast', 'FontSize', 14) 
            grid on 
            hold off
        h = subplot(2,2,3)
            it = 2;
            semilogx(S{it},10^6*K_eff{it},'-o', 'LineWidth',6-it, 'Color', col(2,:))
            hold on
            grid on
            while(it < num_results)
                it = it + 1;
                semilogx(S{it},10^6*K_eff{it},'-o', 'LineWidth',6-it, 'Color', col(it,:))
            end
            xlabel('Nutrient concentration [mM]')
            xticks([10^(-6) 10^(-4) 10^(-2) 10^(0) 10^(2)])
            ylabel('K_{eff} [nM] ', 'FontSize', 16)
            %yticks([0 10 20 30 40 50 60 70 80 90 100 200])
            %ylim([0 200])
            title('Optimal half-saturation')
            set(findall(gca,'-property','FontSize'),'FontSize',16)
            legend(labels(2:end), 'Location', 'southeast', 'FontSize', 14)



            hold off

            subp = subplot(2,2,4);
            it = 2;
            semilogx(S{it}, vol_ratio_periplasm{it},'-o', 'LineWidth',6-it, 'Color', col(2,:))
            hold on
            grid on
            while(it < num_results)
                it = it + 1;
                semilogx(S{it}, vol_ratio_periplasm{it},'-o', 'LineWidth',6-it, 'Color', col(it,:))
            end
            ylabel('Fraction of total volume')
            xlabel('Nutrient concentration [mM]')
            %yticks([0 2 4 6 8 10])
            xticks([10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
            title('Optimal periplasmic volume fraction')
            set(findall(gca,'-property','FontSize'),'FontSize',16)
            legend(labels(2:end), 'Location', 'northeast', 'FontSize', 14)
            hold off

            set(gcf,'Position', [10 10 902 749])
            savefig(strcat(fig_filename, '_Fig1'))

            it = 3
            figure()
            while(it <= length(results)) 

                    subp = subplot(2,2,2*it-5)
                    ax = subp
                    %loglog(S{1}, phi_T{1}+phi_BP{1},'-','LineWidth',8, 'Color',[.98 .6 .51])
                    loglog(S{1}, phi_T{2},'-','LineWidth',4,'Color',col(2,:))
                    hold on
                    %loglog(S{it}, phi_T{it}+phi_BP{it},'-','LineWidth',8, 'Color',[.98 .6 .51])
                    %loglog(S{1}, phi_T{1},'-','LineWidth',4,'Color',col(1,:))
                    loglog(S{it}, phi_T{it},'-','LineWidth',4, 'Color',col(it,:))
                    loglog(S{1}, phi_BP{2}, '-.','LineWidth',4, 'Color', col(2,:))
                    loglog(S{it}, phi_BP{it},'-.','LineWidth',4, 'Color',col(it,:))
                    loglog(S{1}, phi_E{2},':','LineWidth',2, 'Color',col(2,:))
                    loglog(S{it}, phi_E{it},':','LineWidth',2, 'Color',col(it,:))
                    loglog(S{1}, phi_R{2}, '--','LineWidth',2, 'Color', col(2,:))
                    loglog(S{it}, phi_R{it}, '--','LineWidth',2, 'Color', col(it,:))
                    xlim([10^-6 10^2])
                    xticks([10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
                    xlabel('Nutrient concentration [mM]', 'FontSize', 16)
                    ylabel('Fraction of proteome','FontSize', 16)
                    title('Optimal proteome allocation', 'FontSize', 16)
                    set(findall(gca,'-property','FontSize'),'FontSize',16)
                    dummyax = axes('parent',gcf,'position',get(ax,'position'));
                    hold(dummyax, 'on');
                    set(dummyax,'visible','off')
                    h1 = plot(dummyax, nan(1,2), 'LineWidth',3, 'Color', col(2,:));
                    h2 = plot(dummyax, nan(1,2), 'LineWidth',3, 'Color', col(it,:));
                    legend([h1 h2],labels{2},labels{it},'Location', 'south', 'FontSize', 14)
                    dummyax2 = axes('parent',gcf,'position',get(ax,'position'));
                    hold(dummyax2, 'on');
                    set(dummyax2,'visible','off')
                    h7 = plot(dummyax2, nan(1,2), '-', 'LineWidth', 8, 'Color', [0.5 0.5 0.5]);
                    h4 = plot(dummyax2, nan(1,2), '-','LineWidth',4, 'Color', [0 0 0]);
                    h5 = plot(dummyax2, nan(1,2), '-.','LineWidth',4, 'Color', [0 0 0]);
                    h6 = plot(dummyax2, nan(1,2), ':','LineWidth',2, 'Color', [0 0 0]); 
                    h8 = plot(dummyax2, nan(1,2), '--', 'LineWidth',2,'Color', [0 0 0]);
                    legend([h4 h5 h6 h8],'Transport unit', 'Binding protein', 'Metabolism', 'Ribosomes', 'Location','southeast', 'FontSize', 12)
                    grid on 
                    hold off
                subplot(2,2,2*it-4)
                     semilogx(S{2}, ratio_num_BP_num_T{2},'-o', 'LineWidth',5, 'Color', col(2,:))
                    hold on
                        semilogx(S{it}, ratio_num_BP_num_T{it},'-o', 'LineWidth',4, 'Color', col(it,:))
                    ylabel('Ratio')
                    xlabel('Nutrient concentration [mM]')
                    title('Ratio of binding proteins to transport units')
                    %yticks([0 2 4 6 8 10])
                    xticks([10^(-6) 10^(-4) 10^(-2) 10^0 10^2])
                    set(findall(gca,'-property','FontSize'),'FontSize',16)
                    hold off
                    legend(labels{2}, labels{it} , 'Location', 'northeast', 'FontSize', 14)



                it = it + 1;

            end
                   set(gcf,'Position', [10 10 902 749])
                   savefig(strcat(fig_filename, '_Fig2'))
end

    
    
    








