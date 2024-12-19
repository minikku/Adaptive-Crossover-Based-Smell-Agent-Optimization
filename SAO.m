% Smell Agent Optimization (SAO)

% CITATION:
% Salawudeen, A. T., Mu’azu, M. B., Sha’aban, Y. A., & Adedokun, A. E. (2021). A Novel Smell Agent Optimization (SAO):
% An extensive CEC study and engineering application. Knowledge-Based Systems, 232(107486), 107486.
% https://doi.org/10.1016/j.knosys.2021.107486

function [best_cost,best_x,convergence_curve,cul_solution,cul_nfe] = SAO(SearchAgents_no,max_nfe,lb,ub,dim,fobj)
%disp('SAO is now tackling your problem');

rng(sum(100*clock));

% Evaluation settings
eval_setting.case = fobj;
eval_setting.lb = lb;
eval_setting.ub = ub;
eval_setting.dim = dim;
cul_solution = [];
cul_nfe = [];
convergence_curve = [];

% Parameters
nVar = dim;
nPop = SearchAgents_no; % Number of Smell Molecules
m = 2.4; % mass of molecules.
K = 1.38064852*10^(-23); % The boltzman constant
SN = 2.5;
T = 3; %Temperature of gas molecules.
olf = 3.5;

% Main Loop

% Molecules Initialization
molecules = Initialize(nPop,nVar,lb,ub);
v = molecules*0.1;
molecules = molecules + v;

for i=1:nPop
    % Fit(i)=objFun(fitness,molecules(i,:));
    feasible_sol = FeasibleFunction(molecules(i,:),eval_setting);
    [Fit(i),~,~] = feval(fobj,feasible_sol);
    molecules(i,:) = feasible_sol;
end

[Fits,index] = min(Fit); % Obtain the fitness of the best molecule
x_agent = molecules(index,:); % Determin the agent

nfe = 0;
itr = 0;

while nfe < max_nfe

    itr = itr + 1;

    for i=1:nPop
        for j=1:nVar
            % Update the molecular Velocity
            v(i,j)=(v(i,j)+rand*sqrt(3*K*T/m));
        end
    end

    % Implementing the sniffing mode
    for i=1:nPop
        for j=1:nVar
            % molecules(i,j) = molecules(i,j)+v(i,j)+rand*sqrt(3*K*T/m);
            molecules(i,j) = molecules(i,j)+v(i,j);
        end
    end

    for i=1:nPop
        % Fitsniff(i)=objFun(fitness,molecules(i,:));
        feasible_sol = FeasibleFunction(molecules(i,:),eval_setting);
        [Fitsniff(i),~,~] = feval(fobj,feasible_sol);
        molecules(i,:) = feasible_sol;
        nfe = nfe + 1;
    end

    [Fitsniffmin,sindex] = min(Fitsniff);
    xs_agent = molecules(sindex,:);
    [Fitsniffmax,sidx] = max(Fitsniff);
    x_worst = molecules(sidx,:);

    if Fitsniffmin < Fits
        xs_agent = molecules(sindex,:);
        Fits = Fitsniffmin;
    end

    % Evaluate the Trailing mode
    for i=1:nPop
        for j=i:nVar
            molecules(i,j)=molecules(i,j)+rand*olf*(x_agent(1,j)-abs(molecules(i,j)))...
                -rand*olf*(x_worst(1,j)-abs(molecules(i,j)));
        end
    end

    % Make sure no smell molecules excape the boundary
    molecules = clipboundry(lb,ub,molecules);

    % Evaluate the fitness of the Trailing mode
    for i=1:nPop
        % Fitstrail(i)=objFun(fitness,molecules(i,:));
        feasible_sol = FeasibleFunction(molecules(i,:),eval_setting);
        [Fitstrail(i),~,~] = feval(fobj,feasible_sol);
        molecules(i,:) = feasible_sol;
        nfe = nfe + 1;
    end

    [Fitstrailmin,tindex] = min(Fitstrail);

    % Compare the fitness of the trailing mode and the sniffing mode
    % and implement the random mode

    for i=1:nPop
        for j=1:nVar
            if Fitstrail(i) < Fitsniff(i)
                Best_Molecule(i,j) = Fitstrail(1,1);
            elseif Fitstrail(i) > Fitsniff(i)
                molecules(i,j) = molecules(i,j)+rand()*SN;
                % Uncomment the following, if you want to repart sniffing and traling search as a random function.
                % olecules(i,j)=molecules(i,j)+(v(i,j)+rand*sqrt(3*K*T/m));
                % molecules(i,j)=molecules(i,j)+rand*olf*(x_agent(1,j)-abs(molecules(i,j)))...
                % -rand*olf*(x_worst(1,j)-abs(molecules(i,j)));
            end
        end
    end

    for i=1:nPop
        % ybest(i)=objFun(fitness,molecules(i,:));
        feasible_sol = FeasibleFunction(molecules(i,:),eval_setting);
        [ybest(i),~,~] = feval(fobj,feasible_sol);
        molecules(i,:) = feasible_sol;
        nfe = nfe + 1;
    end

    %[SmellObject,Position] = min(ybest);
    %Object(k)=sort(SmellObject,'descend');

    [bestFit,bestIdx] = min(ybest);
    
    best_x = molecules(bestIdx,:);
    best_cost = bestFit;
    convergence_curve = [convergence_curve best_cost];
    cul_solution = [cul_solution best_x];
    cul_nfe = [cul_nfe nfe];
end
end

function mole=Initialize(nPop,nVar,lb,ub)
for i=1:nPop
    for j=1:nVar
        mole(i,j)=lb(j)+rand()*(ub(j)-lb(j));%Defind the initial positions of the smell molecules.
    end
end
end

function molecules=clipboundry(lb,ub,molecules)
[nPop,nVar]=size(molecules);

for i=1:nPop
    for j=1:nVar
        if molecules(i,j)<lb(j)
            molecules(i,j)=lb(j);
        elseif molecules(i,j)>ub(j)
            molecules(i,j)=ub(j);
        end
    end
end
end