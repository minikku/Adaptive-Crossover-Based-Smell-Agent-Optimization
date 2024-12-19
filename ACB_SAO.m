
% =============================================================================================================================================================================
% Adaptive Crossover-Based Smell Agent Optimization (ACB-SAO)

% CITATION:
% P. Duankhan, K. Sunat and C. Soomlek, "An Adaptive Smell Agent Optimization with Binomial Crossover and Linnik Flight for Engineering Optimization Problems," 2024 28th International Computer Science and Engineering Conference (ICSEC), Khon Kaen, Thailand, 2024, pp. 1-6, doi: 10.1109/ICSEC62781.2024.10770710.% =============================================================================================================================================================================

function [best_cost,best_x,convergence_curve,cul_solution,cul_nfe] = ACB_SAO(SearchAgents_no,max_nfe,lb,ub,dim,fobj)

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
% itr = 0;

% Initial CR
CR = 0.9 * ones(nPop, 1);
tau2 = 0.1;
D = dim;

% Golden ratio
golden_ratio = 2/(1 + sqrt(5));

while nfe < max_nfe

    % itr = itr + 1;

    CRold = CR;

    % Sniffing mode
    for i=1:nPop

        % Generate CR
        if rand < tau2
            CR(i) = 0.0 + 1.0 * rand(1, 1);
        end

        trial(i,:) = molecules(i,:);

        % binomial crossover
        jrand = randi(D);
        for j=1:nVar
            if rand <= CR(i) || j == jrand
                trial(i,j) = molecules(i,j)+v(i,j)+rand*sqrt(3*K*T/m);
                % trial(i,j) = molecules(i,j)+v(i,j);
            end
        end

        % % exponential crossover
        % jrand = randi(D);
        % L = 0;
        % while (rand <= CR(i)) && (L < D)
        %     j = mod(jrand + L - 1, D) + 1;
        %     trial(i,j) = molecules(i,j)+v(i,j)+rand*sqrt(3*K*T/m);
        %     L = L + 1;
        % end

        % Make sure no smell molecules excape the boundary
        trial(i,:) = clipboundry(lb,ub,trial(i,:));

        % Evaluation
        feasible_sol = FeasibleFunction(trial(i,:),eval_setting);
        [Fitsniff(i),~,~] = feval(fobj,feasible_sol);
        trial(i,:) = feasible_sol;
        nfe = nfe + 1;

        if Fitsniff(i) <= Fit(i)
            molecules(i,:) = trial(i,:);
            Fit(i) = Fitsniff(i);
        else
            CR(i) = CRold(i);
        end
    end

    % for i=1:nPop
    %     for j=1:nVar
    %         % Update the molecular Velocity
    %         v(i,j)=(v(i,j)+rand*sqrt(3*K*T/m));
    %     end
    % end

    [Fitsniffmin,sindex] = min(Fitsniff);
    xs_agent = molecules(sindex,:);
    [~,sidx] = max(Fitsniff);
    x_worst = molecules(sidx,:);

    CRold = CR;

    % Trailing mode
    for i=1:nPop

        % Generate CR
        if rand < tau2
            CR(i) = 0.0 + 1.0 * rand(1, 1);
        end

        trial(i,:) = molecules(i,:);

        % for j=i:nVar
        %     % molecules(i,j)=molecules(i,j)+rand*olf*(x_agent(1,j)-abs(molecules(i,j)))...
        %     %     -rand*olf*(x_worst(1,j)-abs(molecules(i,j)));
        %     mutant(i,j)=molecules(i,j)+rand*olf*(x_agent(1,j)-abs(molecules(i,j)))...
        %         -rand*olf*(x_worst(1,j)-abs(molecules(i,j)));
        % end

        % binomial crossover
        jrand = randi(D);
        for j=1:nVar
            if rand <= CR(i) || j == jrand
                trial(i,j) = molecules(i,j)+rand*olf*(x_agent(1,j)-abs(molecules(i,j)))-rand*olf*(x_worst(1,j)-abs(molecules(i,j)));
            end
        end

        % % exponential crossover
        % jrand = randi(D);
        % L = 0;
        % while (rand <= CR(i)) && (L < D)
        %     j = mod(jrand + L - 1, D) + 1;
        %     trial(i,j) = molecules(i,j)+rand*olf*(x_agent(1,j)-abs(molecules(i,j)))-rand*olf*(x_worst(1,j)-abs(molecules(i,j)));
        %     L = L + 1;
        % end

        % Make sure no smell molecules excape the boundary
        trial(i,:) = clipboundry(lb,ub,trial(i,:));

        % Evaluation
        feasible_sol = FeasibleFunction(trial(i,:),eval_setting);
        [Fitstrail(i),~,~] = feval(fobj,feasible_sol);
        trial(i,:) = feasible_sol;
        nfe = nfe + 1;

        if Fitstrail(i) <= Fit(i)
            molecules(i,:) = trial(i,:);
            Fit(i) = Fitstrail(i);
        else
            CR(i) = CRold(i);
        end
    end

    CRold = CR;

    % Random mode
    for i=1:nPop

        % Generate CR
        if rand < tau2
            CR(i) = 0.0 + 1.0 * rand(1, 1);
        end

        trial(i,:) = molecules(i,:);

        % binomial crossover
        jrand = randi(D);
        for j=1:nVar
            if rand <= CR(i) || j == jrand
                trial(i,j) = molecules(i,j)+LnF2(golden_ratio,0.05,1,1);
            end
        end

        % % exponential crossover
        % jrand = randi(D);
        % L = 0;
        % while (rand <= CR(i)) && (L < D)
        %     j = mod(jrand + L - 1, D) + 1;
        %     trial(i,j) = molecules(i,j)+LnF2(golden_ratio,0.05,1,1);
        %     L = L + 1;
        % end

        % Make sure no smell molecules excape the boundary
        trial(i,:) = clipboundry(lb,ub,trial(i,:));

        % Evaluation
        feasible_sol = FeasibleFunction(trial(i,:),eval_setting);
        [Fitsrand(i),~,~] = feval(fobj,feasible_sol);
        trial(i,:) = feasible_sol;
        nfe = nfe + 1;

        if Fitsrand(i) <= Fit(i)
            molecules(i,:) = trial(i,:);
            Fit(i) = Fitsrand(i);
        else
            CR(i) = CRold(i);
        end
    end

    [bestFit,bestIdx] = min(Fit);

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
        mole(i,j)=lb(j)+rand()*(ub(j)-lb(j)); % Defind the initial positions of the smell molecules.
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

function trial = binomialCrossover(target, mutant, CR, D)
jrand = randi(D);
trial = target;
for j = 1:D
    if rand <= CR || j == jrand
        trial(j) = mutant(j);
    end
end
end

function trial = exponentialCrossover(target, mutant, CR, D)
jrand = randi(D);
L = 0;
trial = target;
while (rand <= CR) && (L < D)
    trial(mod(jrand + L - 1, D) + 1) = mutant(mod(jrand + L - 1, D) + 1);
    L = L + 1;
end
end

function xhold = LnF2(alpha,scale,m, n)
xhold = laprnd(m, n, 0, 1);
% xhold = randl(m,n); %Laplacian distributed pseudorandom numbers.
SE = sign(rand(m,n)-0.5) .* xhold;
% U = rand(1,n);
U = rand(m,n);
xhold = (sin(0.5*pi*alpha).*cot(0.5*pi*(alpha*U))-cos(0.5*pi*alpha)).^(1/alpha);
xhold = scale .* SE ./ xhold;
end

function y = laprnd(m, n, mu, sigma)
% LAPRND generate i.i.d. laplacian random number drawn from laplacian distribution
%   with mean mu and standard deviation sigma.
%   mu      : mean
%   sigma   : standard deviation
%   [m, n]  : the dimension of y.
%   Default mu = 0, sigma = 1.
%   For more information, refer to
%   http://en.wikipedia.org./wiki/Laplace_distribution

%   Author  : Elvis Chen (bee33@sjtu.edu.cn)
%   Date    : 01/19/07

% Check inputs
if nargin < 2
    error('At least two inputs are required');
end

if nargin == 2
    mu = 0; sigma = 1;
end

if nargin == 3
    sigma = 1;
end

% Generate Laplacian noise
u = rand(m, n)-0.5;
b = sigma / sqrt(2);
y = mu - b * sign(u).* log(1- 2* abs(u));
end