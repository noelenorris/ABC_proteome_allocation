%% Km plot

clear all
%close all

% colormap
col = lines(7);

%ABC kinetics rates
k0r = 100/1000; %msec-1
k0f = 10^5/1000; %mM-1msec-1
k1f = 0.01*21; %mM-1msec-1
k2 = 210/1000; %msec-1
k3 = (210/1000)*0.01; %msec-1




        BP_total_it = 10.^[-4:.1:2];
        T_total_it = [.01 .1 1.16];
                
       for i=1:1:length(T_total_it)
           for j=1:1:length(BP_total_it)
                   
                   T_total = T_total_it(i);
                   BP_total = BP_total_it(j);


                    %Calculuate V_max 
                    S_p2 = 10^10;


                    T_S_BP = (k3.*(k2.*k3.*k0r - (BP_total.^2.*S_p2.^2.*k2.^2.*k0f.^2.*k1f.^2 + 2.*BP_total.^2.*S_p2.^2.*k2.*k3.*k0f.^2.*k1f.^2 + BP_total.^2.*S_p2.^2.*k3^2*k0f^2*k1f^2 ...
                           - 2*BP_total.*S_p2.^2.*T_total.*k2^2.*k0f^2.*k1f^2 - 4.*BP_total.*S_p2^2.*T_total*k2*k3*k0f^2*k1f^2 - 2*BP_total.*S_p2.^2.*T_total.*k3^2*k0f^2*k1f^2 ...
                           + 2*BP_total.*S_p2.^2.*k2^2*k3*k0f^2*k1f + 2*BP_total.*S_p2.^2.*k2*k3^2*k0f^2*k1f - 2*BP_total.*S_p2.*T_total.*k2^2*k3*k0f*k1f^2 - 2*BP_total.*S_p2.*T_total.*k2*k3^2*k0f*k1f^2 ...
                           + 2*BP_total.*S_p2.*k2^2*k3*k0f*k1f*k0r + 2*BP_total.*S_p2*k2*k3^2*k0f*k1f*k0r + S_p2.^2.*T_total.^2*k2^2*k0f^2*k1f^2 + 2*S_p2.^2.*T_total.^2*k2*k3*k0f^2*k1f^2 ...
                           + S_p2.^2.*T_total.^2*k3^2*k0f^2*k1f^2 + 2.*S_p2.^2.*T_total*k2^2*k3*k0f^2*k1f + 2*S_p2^2*T_total*k2*k3^2*k0f^2*k1f + S_p2^2*k2^2*k3^2*k0f^2 + 2*S_p2.*T_total.^2*k2^2*k3*k0f*k1f^2 ...
                           + 2*S_p2.*T_total.^2*k2*k3^2*k0f*k1f^2 + 2*S_p2.*T_total*k2^2*k3^2*k0f*k1f + 2*S_p2.*T_total*k2^2*k3*k0f*k1f*k0r + 2*S_p2.*T_total*k2*k3^2*k0f*k1f*k0r + 2*S_p2.*k2^2*k3^2*k0f*k0r ...
                           + T_total.^2*k2^2*k3^2*k1f^2 + 2*T_total.*k2^2*k3^2*k1f*k0r + k2^2*k3^2*k0r^2).^(1/2) + S_p2.*k2*k3*k0f + T_total.*k2*k3*k1f + BP_total.*S_p2*k2*k0f*k1f + BP_total.*S_p2*k3*k0f*k1f ...
                           + S_p2.*T_total*k2*k0f*k1f + S_p2.*T_total*k3*k0f*k1f))./(2*k1f*(k2 + k3)*(k2*k3 + S_p2.*k2*k0f + S_p2.*k3*k0f));


                     % cytoplasmic transport rate ([S]_p --> S_c)
                     Vmax = k2*T_S_BP;


                   S_p2 = 10.^[-10:.01:-2];
                   
                    T_S_BP = (k3.*(k2.*k3.*k0r - (BP_total.^2.*S_p2.^2.*k2.^2.*k0f.^2.*k1f.^2 + 2.*BP_total.^2.*S_p2.^2.*k2.*k3.*k0f.^2.*k1f.^2 + BP_total.^2.*S_p2.^2.*k3^2*k0f^2*k1f^2 ...
                   - 2*BP_total.*S_p2.^2.*T_total.*k2^2.*k0f^2.*k1f^2 - 4.*BP_total.*S_p2.^2.*T_total*k2*k3*k0f^2*k1f^2 - 2*BP_total.*S_p2.^2.*T_total.*k3^2*k0f^2*k1f^2 ...
                   + 2*BP_total.*S_p2.^2.*k2^2*k3*k0f^2*k1f + 2*BP_total.*S_p2.^2.*k2*k3^2*k0f^2*k1f - 2*BP_total.*S_p2.*T_total.*k2^2*k3*k0f*k1f^2 - 2*BP_total.*S_p2.*T_total.*k2*k3^2*k0f*k1f^2 ...
                   + 2*BP_total.*S_p2.*k2^2*k3*k0f*k1f*k0r + 2*BP_total.*S_p2*k2*k3^2*k0f*k1f*k0r + S_p2.^2.*T_total.^2*k2^2*k0f^2*k1f^2 + 2*S_p2.^2.*T_total.^2*k2*k3*k0f^2*k1f^2 ...
                   + S_p2.^2.*T_total.^2*k3^2*k0f^2*k1f^2 + 2.*S_p2.^2.*T_total*k2^2*k3*k0f^2*k1f + 2*S_p2.^2.*T_total*k2*k3^2*k0f^2*k1f + S_p2.^2*k2^2*k3^2*k0f^2 + 2*S_p2.*T_total.^2*k2^2*k3*k0f*k1f^2 ...
                   + 2*S_p2.*T_total.^2*k2*k3^2*k0f*k1f^2 + 2*S_p2.*T_total*k2^2*k3^2*k0f*k1f + 2*S_p2.*T_total*k2^2*k3*k0f*k1f*k0r + 2*S_p2.*T_total*k2*k3^2*k0f*k1f*k0r + 2*S_p2.*k2^2*k3^2*k0f*k0r ...
                   + T_total.^2*k2^2*k3^2*k1f^2 + 2*T_total.*k2^2*k3^2*k1f*k0r + k2^2*k3^2*k0r^2).^(1/2) + S_p2.*k2*k3*k0f + T_total.*k2*k3*k1f + BP_total.*S_p2*k2*k0f*k1f + BP_total.*S_p2*k3*k0f*k1f ...
                   + S_p2.*T_total*k2*k0f*k1f + S_p2.*T_total*k3*k0f*k1f))./(2*k1f*(k2 + k3)*(k2*k3 + S_p2.*k2*k0f + S_p2.*k3*k0f));

                    uptake_K_eff = k2*T_S_BP;

                    if(~isnan(uptake_K_eff(1)))

                        f_uptake = fitnlm(S_p2', (uptake_K_eff/Vmax)', 'y ~ x/(b1+x)', 10^(-5));

                        K_M_exact{i,j} = f_uptake.Coefficients.Estimate;
                    else
                        K_M_exact{i,j} = NaN;
                    end
                    
                    K_T = k3*k2/(k1f*(k2+k3));
                    K_D = k0r/k0f;
                    K_approx{j} = K_T*K_D./(K_T + BP_total);
                end
           end
   
       
       
       K_M = cell2mat(K_M_exact);
       K_approx = cell2mat(K_approx);
       
       figure()
        loglog(BP_total_it,10^6*K_approx,':', 'LineWidth',10, 'Color', col(3,:))
        hold on
        loglog(BP_total_it,10^6*squeeze(K_M(1,:)),'-', 'LineWidth',8, 'Color', col(4,:))
        loglog(BP_total_it,10^6*squeeze(K_M(2,:)),'-', 'LineWidth',6, 'Color', col(5,:))
        loglog(BP_total_it,10^6*squeeze(K_M(3,:)),'-', 'LineWidth',4, 'Color', col(6,:))
        plot(7.7, 2.91, 'r*', 'MarkerSize', 15, 'LineWidth', 6)
        xlabel('[BP]_{total} [mM]', 'FontSize', 16)
        %xticks([10^(-4) 10^(-2) 10^(0) 10^(2) 10^(4)])
        ylabel('K_M [nM] ', 'FontSize', 16)
        legend('Approximation', '[T]_{total} = 10 \muM', '[T]_{total} = 100 \muM', '[T]_{total} = 1.16 mM', 'Optimal solution at 1nM')
        %yticks([0 5 10 15 20 25])
        set(findall(gcf,'-property','FontSize'),'FontSize',18)
        grid on
        hold off