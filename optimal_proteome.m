%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                          %
%  OPTIMAL_PROTEOME.M: A proteome allocation problem with nutrient         %
%  transport via an ABC transport system or PTS.                           %
%                                                                          %
%  !IMPORTANT! PLEASE READ: This code will not run  without                %
%       MATLAB's Global Optimization Toolbox. To install the toolbox,      %
%       click: Home --> Add-Ons --> Get Add-Ons.                           %
%       Search for "Global Optimization Toolbox" and click on toolbox      %
%       then "Install". Follow the prompts to finish installation.         %
%       (In older versions of MATLAB, you can obtain the toolbox           %
%        by re-running the installer.)                                     %
%                                                                          %
%   This code iterates over specified optimization problems, calling       %
%   fmincon with constraints specified in helper function                  %
%   eqcon_optimal_proteome.                                                %
%                                                                          %
%   Because this is a nonlinear optimization problem with complex          %
%   parameter landscape, this code solves the same optimization problem    %
%   num_it times, each time with a different random initial condition.     %
%                                                                          %
%      Noele Norris                                                        %
%                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [optimal_values] = optimal_proteome(parameter_mat_name)

    warning('off','MATLAB:nearlySingularMatrix')
     
    file_name = strcat('results_', datestr(now,'mmmmddyyyyHHMMSSFFF'), '_', parameter_mat_name)
    
    %load values specifying protoeme allocation problem
    load(parameter_mat_name)
    
    num = length(x_axis_values);
    optimal_values = cell(num,1);
    
    solutions = [];

    % iterate over different optimization problems
    for p=[1:1:num]
        
        optimal_values{p} = zeros(13,1)';
        
        x_axis_value = x_axis_values(p)

        % specify fmincon optimization problem
        A_op = [];
        b_op = [];
        Aeq_op = [];
        beq_op = [];
       
        % bounds
        lb = min_bounds./z(:,p);
        ub = max_bounds./z(:,p);

        % optimization function, maximize growth rate
        fun = @(x)-log(x(1));

        % helper function specifies constraints of nonlinear optimization problem
        nonlcon = @(x)eqcon_optimal_proteome(x,p,parameter_mat_name);

        % options for optimization
        options = optimoptions('fmincon','Algorithm','interior-point',...
                               'Display','final-detailed',...
                               'OptimalityTolerance',1e-2, ...
                               'ConstraintTolerance',1e-5, ...
                               'MaxFunctionEvaluations',2e4, ...
                               'MaxIterations', 1e4, ...
                               'TolFun', 1e-6);


        x = zeros(13,1);

        tic

        % solve optimization problems with different, random initial conditions
        for i=1:1:num_it

            i  

            x0 = [1*rand(13,1)]

            [x_new,fval,exitflag,output] = fmincon(fun,x0,A_op,b_op,Aeq_op,beq_op,lb,ub,nonlcon, options);

            exitflag;
            output;

            if(exitflag > 0) %if solution is found
                
                s = [x_axis_value; x_new];
                solutions = [solutions s];
                
               if(x_new(1) > x(1)) %if new growth rate is higher, found new optimal solution
                   
                    x_old = x;
                    x = x_new;
                    disp(parameter_mat_name)
                    x_axis_value
                    disp('Current global optimal vs. previous global:')
                    [x x_old]
                else
                    disp(parameter_mat_name)
                    x_axis_value
                    disp('Local optimal vs. global:')
                   [x_new x]
               end
               
            end
            
            toc
        end

        x
        optimal_values{p} = x';

        save(file_name)
    end
    
end