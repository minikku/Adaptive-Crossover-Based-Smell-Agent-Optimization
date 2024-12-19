
% =============================================================================================================================================================================
% Adaptive Crossover-Based Smell Agent Optimization (ACB-SAO)

% CITATION:
% P. Duankhan, K. Sunat and C. Soomlek, "An Adaptive Smell Agent Optimization with Binomial Crossover and Linnik Flight for Engineering Optimization Problems," 2024 28th International Computer Science and Engineering Conference (ICSEC), Khon Kaen, Thailand, 2024, pp. 1-6, doi: 10.1109/ICSEC62781.2024.10770710.% =============================================================================================================================================================================


clc;
clearvars;

all_algo = {
    'ACB_SAO', ...  	1. Adaptive Crossover-Based Smell Agent Optimization (ACB-SAO)
    'SAO', ...  	2. Smell Agent Optimization (SAO)
    };

run_test = 1; % Run test
scoring = 0; % Show scores

selected_test_algo = [1,2];

selected_problem = 6; %[1,3:30];

% Brenchmark
brenchmark_name = 'CEC2017';
all_dim_size = 100; %[10,30,50,100];

if run_test

    % for SearchAgents_no = pop_size % Population size
    %     for Max_iteration = max_nfes % fitness evaluation

    % Brenchmark
    all_func_no = selected_problem; % F1-30

    % Create output directories for the experimental.
    output_dir_name = strcat(brenchmark_name,'/test_');
    OutputDir = strcat('./results/',output_dir_name);

    if ~exist(OutputDir, 'dir')
        mkdir(OutputDir);
    end

    for algo_no = selected_test_algo

        % algorithm's name
        algo_name = string(all_algo(algo_no));

        % Create output directories for the experimental.
        algo_dir_name = strcat(OutputDir,'/',algo_name);

        if ~exist(algo_dir_name, 'dir')
            mkdir(algo_dir_name);
        end

        for dim_size = all_dim_size

            % Create output directories for the experimental.
            dim_dir_name = strcat(algo_dir_name,'/D',num2str(dim_size));

            if ~exist(dim_dir_name, 'dir')
                mkdir(dim_dir_name);
            end

            for func_no = all_func_no

                % Create output directories for the experimental.
                func_no_dir_name = strcat(dim_dir_name,'/F',num2str(func_no));

                if ~exist(func_no_dir_name, 'dir')
                    mkdir(func_no_dir_name);
                end

                % CEC17 F2 has been restored.
                Function_name = strcat('F',num2str(func_no));
                desired_dim = dim_size;

                % Load details of the selected benchmark function
                switch brenchmark_name
                    case 'CEC2005'
                        [lb,ub,dim,fobj]=Get_Functions_details(Function_name,desired_dim);
                    case 'CEC2013'
                        [lb,ub,dim,fobj]=CEC2013_Function(Function_name,desired_dim);
                    case 'CEC2014'
                        [lb,ub,dim,fobj]=CEC2014_Function(Function_name,desired_dim);
                    case 'CEC2017'
                        [lb,ub,dim,fobj]=CEC2017_Function(Function_name,desired_dim);
                end

                lb = lb .* ones(1,dim); % Lower bound
                ub = ub .* ones(1,dim); % Upper bound

                SearchAgents_no = 30;
                Max_iteration = 10000 * dim;

                % Cost Function
                funj = @(x) CostFunction(x, fobj);

                if desired_dim ~= dim
                    fprintf('The operating dimension is set to be %d\n',dim);
                end

                watchRunAll = tic; % Elapsed time for all run.

                parfor run = 1:52

                    watchRun = tic; % Elapsed time for each run.

                    try
                        [best_score(run,:),best_pos(run,:),cg_curve(run,:),~,~] = feval(algo_name,SearchAgents_no,Max_iteration,lb,ub,dim,funj);
                    catch
                        [best_score(run,:),best_pos(run,:),cg_curve(run,:)] = feval(algo_name,SearchAgents_no,Max_iteration,lb,ub,dim,funj);
                    end


                    elapsedRun = toc(watchRun);
                    elapsed_time_run(run,:) = elapsedRun; % Elapsed time for each run.
                end

                fprintf(strcat('✔\t',brenchmark_name,':\t',algo_name,' F',num2str(func_no),' D',num2str(dim_size),'\n'));
                fprintf('Average best score: %e \nMedian: %e \nStd: %e\n\n',mean(best_score(:,1)),median(best_score(:,1)),std(best_score(:,1)));

                elapsedRunAll = toc(watchRunAll); % Elapsed time for all run.

                clear best_score best_pos cg_curve elapsed_time_run cul_solution cul_nfe record_data;
            end
        end

        fprintf(strcat('\n'));

    end

    %     end % max_itr loop
    % end %agent_no loop
end

disp('----- Done -----');
