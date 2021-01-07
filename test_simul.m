clc;clear;close all;
rng('default');
% Simulating dataset
nlv = 3; % Latent variable
x = 1:100; % Spectral wavenumber channel index
T = 10:10:70; % Temperature from 10 to 70 ¡æ
% Concentration profiles
c0 = linspace(0.1,1,20)';
A=[c0.^2 c0.^0.5 c0];

% Spectral profiles whith peak width and  position shift against temperature
B = nan([length(x) nlv length(T)]);
B(:,:,1) = [gaussmf(x,[3,20]);gaussmf(x,[5,45]);gaussmf(x,[10,80])]';
B(:,:,2) = [gaussmf(x,[3,23]);gaussmf(x,[5,47]);gaussmf(x,[10,81])]';
B(:,:,3) = [gaussmf(x,[3,26]);gaussmf(x,[5,49]);gaussmf(x,[10,82])]';
B(:,:,4) = [gaussmf(x,[3,29]);gaussmf(x,[5,51]);gaussmf(x,[10,83])]';
B(:,:,5) = [gaussmf(x,[3,32]);gaussmf(x,[5,53]);gaussmf(x,[10,84])]';
B(:,:,6) = [gaussmf(x,[3,35]);gaussmf(x,[5,55]);gaussmf(x,[10,85])]';
B(:,:,7) = [gaussmf(x,[3,37]);gaussmf(x,[5,57]);gaussmf(x,[10,86])]';

% Temperature profiles
C=[T' 80-T' max([T' 80-T']+10,[],2)];

% Generating the simulated tree-way dataset
X(:,:,1) = A*diag(C(1,:))*squeeze(B(:,:,1))';
X(:,:,2) = A*diag(C(2,:))*squeeze(B(:,:,2))';
X(:,:,3) = A*diag(C(3,:))*squeeze(B(:,:,3))';
X(:,:,4) = A*diag(C(4,:))*squeeze(B(:,:,4))';
X(:,:,5) = A*diag(C(5,:))*squeeze(B(:,:,5))';
X(:,:,6) = A*diag(C(6,:))*squeeze(B(:,:,6))';
X(:,:,7) = A*diag(C(7,:))*squeeze(B(:,:,7))';

% Adding random noise
noise = 0.001*max(X(:))*rand(size(X));
X = X+noise;
DimX=size(X);
clear c0 noise T x
%% Non-negative constraint
dntd_non_negative_const=[2 2 2];% Non-negative constraint for three profiles
[A_ini,B_ini,C_ini] = dtld(X,nlv);
% A_ini=rand(DimX(1),nlv);
% B_ini=rand(DimX(2),nlv);
% C_ini=rand(DimX(3),nlv);
[Factor_dntd] = dntd(X,nlv,dntd_non_negative_const,'Thres',1e-6,'Niter',200,'NStepNoLearning',0,'SmoothB',0,'Init',{A_ini,B_ini,C_ini});
[A_dntd,B_dntd,C_dntd]= deal(Factor_dntd{1},Factor_dntd{2},Factor_dntd{3});

%ploting the simulated profiles
plot3way(A,B,C);
%plot resolved profiles by the proposed method
plot3way(A_dntd,B_dntd,C_dntd);