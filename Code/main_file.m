% Each section can be run independtly

%% Import & Data cleaning
clear; close all; clc;

% Import investment data in a 1xT array. Each cell is a yearly NxN matrix
investment_network=sheet_import('INVEST_network.xlsx');

% Same for IO flows
prod_flows=sheet_import('PROD_flows.xlsx');

% Import sectors data (TxN)
sectors_data=xlsread('sectors_data_quarterly.xlsx');

% Need to convert the flows to a weights: x_ijt / SUM_j x_ijt (input divided by 
% total input expenditure)
for t = 1:72
    A=prod_flows{t};
    col_sum=sum(A);
    for j = 1:37
        for i = 1:37
            A(i,j)=A(i,j)/col_sum(j);
        end
    end
    production_network{t}=A;
end

% Import money supply data
MZM=readmatrix("MZM_data.csv");
Msize=size(MZM);

% Make money supply quarterly (end of the 3 months) 
index=1;
for q=1:Msize(1,2)/3
    MZM3(1,q)=MZM(1,index+2);
    index=index+3;
end

% Transform money supply into money supply shocks (M_t - M_{t-1} / M_{t-1})
for q=2:Msize(1,2)/3
    M(1,q-1)=(MZM3(q)-MZM3(q-1))/MZM3(q-1);
end

% Keep only the Q1 2005 - Q4 2018 data
sectors_data=sectors_data(1:60,:);
M=M(1,184:243);

% Export data
clear i k sheet_name A t i j col_sum Msize MZM MZM3 q index;
save('data.mat')

%% FIGURE 1: US IO graph 2012

clear; close all; clc;
data=csvread('usinputoutput2012_200.csv',1,1);

% Custom colormap if wanted
%colormapeditor

% Create a matrix of start points (S), end points (E), and weight for edges
edgeC = [.7 .7 .7];
nodeC = [0 0 .26];

% Plot the graph
G=digraph(data);
p=plot(G,'NodeColor',nodeC,'EdgeAlpha',0.13,'MarkerSize',3,'LineWidth',1);

%% FIGURE 3: US IO+INVEST graph 2018

clear; close all; clc;

data=csvread('Investmentnetwork_2018.csv');
data2=csvread('Productionnetwork_2018.csv');

% Custom colormap if wanted
%colormapeditor

edgeBLUE = [.7 .7 .7];
nodeBLUE = [0 0 .26];

edgeRED = [1 0 0];
nodeRED = [0.6350 0.0780 0.1840];

% 2018 investment Network
G=digraph(data, "omitselfloops");
p=plot(G,'NodeColor',nodeRED, 'EdgeColor',edgeRED,'LineWidth',0.9,'EdgeAlpha',0.2,'MarkerSize',3,'Layout','force');

% 2018 production network
G=digraph(data2, "omitselfloops");
p=plot(G,'NodeColor',nodeBLUE, 'EdgeAlpha',0.2,'MarkerSize',3,'LineWidth',0.9,'Layout','force');

%% TABLE 1: Stationarity of IO matrices
clear; close all; clc;
load("data.mat");

% Parameters (to make less computing-intensive tests)
N=37; % Max: 37
T=72; % Max: 72
r=59; % Reduced sample, r=59 starts at 2005. NB: we must have T>r

% Vectorize the networks so that each row is a coefficient time series
A_ts=vectorize_IO(production_network,N,T);
Theta_ts=vectorize_IO(investment_network,N,T);

% ADF test of unit root
[ADFA.not, ADFA.ten, ADFA.five, ADFA.one]=matrix_ADFtest(A_ts,N,T); % A, full sample
[ADFAr.not, ADFAr.ten, ADFAr.five, ADFAr.one]=matrix_ADFtest(A_ts(:,r:T),N,T-r); % A, 2005-2018
[ADFT.not, ADFT.ten, ADFT.five, ADFT.one]=matrix_ADFtest(Theta_ts,N,T); % Theta, full sample
[ADFTr.not, ADFTr.ten, ADFTr.five, ADFTr.one]=matrix_ADFtest(Theta_ts(:,r:T),N,T-r); % Theta, 2005-2018

%% TABLE 2: regression 1
clear; close all; clc;
load("data.mat");

% Parameters (to make less computing-intensive tests)
N=36; % Max: 37 => need to remove one otherwise we have perfect colinearity
T=60; % Max: 60 quarters
ref=65; % Reference year for IO matrices. 59=2005, 65=2011 72=2018

% Matrix Y: TxN gross output over time
for t=1:T
    for i=1:N
        Y(t,i)=sectors_data(t,i);
    end
end

% Reference IO linkages
A=production_network{ref}(1:N,1:N);% Input-output matrix (materials)
Lm=inv(eye(N)-A); % Leontief inverse (materials)
Theta=investment_network{ref}(1:N,1:N);% Input-output matrix (capital)
Lk=inv(eye(N)-Theta); % Leontief inverse (capital)

% Matrix X: TxN
for t=1:T
    xi=Lm.*M(t); % Leontief * shock
    X{t}= [eye(N),xi]; % Add a constant for each sector (sectoral FE)
end

X=transpose(X);

% Regression
[beta,Sigma]=mvregress(X,Y,'maxiter',10000);
alpha=beta(1:36,1);
beta=beta(37:72,1); % remove alpha

std_error=sqrt(diag(Sigma));
direct_effect=beta;
total_effect=Lm*direct_effect;
indirect_effect=total_effect-direct_effect;
network_share=indirect_effect./total_effect;


%% TABLE 2: regression 2
clear; close all; clc;
load("data.mat");

% Parameters (to make less computing-intensive tests)
N=36; % Max: 37 => need to remove one otherwise we have perfect colinearity
T=60; % Max: 60 quarters
ref=65; % Reference year for IO matrices. 59=2005, 65=2011 72=2018

% Matrix Y: TxN gross output over time
for t=1:T
    for i=1:N
        Y(t,i)=sectors_data(t,i);
    end
end

% Reference IO linkages
A=production_network{ref}(1:N,1:N);% Input-output matrix (materials)
Lm=inv(eye(N)-A); % Leontief inverse (materials)
Theta=investment_network{ref}(1:N,1:N);% Input-output matrix (capital)

% Populate the 0 values of Capital Leontiev with random errors to avoid a
% rank defficient matrix (consistent with Foerster 2011)
for i=1:N
    for j=1:N
        if Theta(i,j)==0
            Theta(i,j)=normrnd(0.0001,0.0001);
        end
    end
end 

Lk=inv(eye(N)-Theta); % Leontief inverse (capital)

% Add the two Leontiev inverse
Ltot=Lm+Lk;

% Matrix X: TxN
for t=1:T
    xi=Ltot.*M(t); % Leontief * shock
    X{t}= [eye(N),xi]; % Add a constant for each sector (sectoral FE)
end

X=transpose(X);

% Regression
[beta,Sigma]=mvregress(X,Y,'maxiter',10000);
alpha=beta(1:36,1);
beta=beta(37:72,1); % remove alpha

std_error2=sqrt(diag(Sigma));
direct_effect2=beta;
total_effect2=Ltot*direct_effect2;
indirect_effect2=total_effect2-direct_effect2;
network_share2=indirect_effect2./total_effect2;


%% TABLE 3: Acemoglu upstream/downstream regression
clear; close all; clc;
load("data.mat");

% Parameters (to make less computing-intensive tests)
N=36; % Max: 37 => need to remove one otherwise we have perfect colinearity
T=60; % Max: 60 quarters, where the first observation is droped for the AR(1)
ref=72; % Reference year for IO matrices. 59=2005, 65=2011 72=2018

% Reference IO linkages
A=production_network{ref}(1:N,1:N);% Input-output matrix (materials)
Lm=inv(eye(N)-A); % Leontief inverse (materials)
Theta=investment_network{ref}(1:N,1:N);% Input-output matrix (capital)

% Populate the 0 values of Capital Leontiev with random errors to avoid a
% rank defficient matrix (consistent with Foerster 2011)
for i=1:N
    for j=1:N
        if Theta(i,j)==0
            Theta(i,j)=normrnd(0.0001,0.0001);
        end
    end
end 
Lk=inv(eye(N)-Theta); % Leontief inverse (capital)

% Matrix Y: TxN gross output over time
for t=2:T
    for i=1:N
        Y(t-1,i)=sectors_data(t,i);
    end
end

% Matrices Y(t-1): TxN gross output over time
for t=1:T-1
    for i=1:N
        Yt1(t,i)=sectors_data(t,i);
    end
end

% Compute the upstream and downstream TxN matrices
Ltot=Lm+Lk; 
for t=1:T
    for i=1:N
        colsum=sum(Ltot,1);
        rowsum=sum(Ltot,2);
        down(t,i)=(colsum(1,i)-Ltot(i,i))*M(1,t);
        up(t,i)=(rowsum(i,1)-Ltot(i,i))*M(1,t);
    end 
end

% Matrix X: 59 times a NxN^4 matrix
for t=2:T
    cst=ones(N,1); % Add a constant for each sector (sectoral FE)
    x1=transpose(Yt1(t-1,:));
    x2=transpose(up(t,:)); % Up shocks
    x3=transpose(down(t,:));
    X{t-1}= [cst,x1,x2,x3]; % concatenate
end
X=transpose(X);

% Regression
[beta,Sigma,E,CovB,logL]=mvregress(X,Y,'maxiter',10000);
omega=beta(1,1);
psi=beta(2,1);
beta_up=beta(3,1);
beta_down=beta(4,1);
std_error3=sqrt(diag(CovB));
