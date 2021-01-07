function [Factor,Results] = dntd(X,F,const,varargin)
% dntd direct non-trilinear decomposition
% ----------------------INPUT---------------------
% X          X is the input array, which can be from three- to N-way (also
%            twoway if the third mode is interpreted as a onedimensional
%            mode).
%
% F          Number of factors/components sought.
% ----------------OPTIONAL INPUT---------------------
% const      A vector telling type of constraints put on the loadings of the
%            different modes. Same size as DimX but the i'th element tells
%            what constraint is on that mode.
%            0 => no constraint,
%            2 => nonnegativity
%            6 => gauss profile constrain (Only applicable for shifing B model)
%            7 => gauss and Lorentzian profile constrain for B model
%            For no constraints in a threeway problem const = [0 0 0]
% Options    Optional parameters. If not given or set to zero or [],
%            defaults will be used. If you want Options(5) to be 2 and
%            not change others, simply write Options(5)=2. Even if Options
%            hasn't been defined Options will contain zeros except its
%            fifth element.
%
%            Thres - Convergence criterion
%            The relative change in fit for which the algorithm stops.
%            Standard is 1e-6, but difficult data might require a lower value.
%
%            Init - Initialization method
%            If the Init is set to a cell with three matrix, it will be set
%            as the initial values. Otherwise, the initional will be set to
%            random values.
%
% Written by J zhang
% zhangjin@mail.nankai.edu.cn

tic;
addpath('./N_way_toolbox')
if rem(nargin-3,2);error('Parameters to dntd should be pairs');end
if nargin<3;const=[2 2 2];end
Niter=200;     % The default number of iteration
Thres=1e-9;    % The default stop condition for iteration
SmoothB = true; % The default choice for profile smoothing
numParameters = (nargin-3)/2;
NStepNoLearning=5; % The default number of steps no imposing average constraint
Init = 'random';
for n =1:numParameters
    Parameters = varargin{(n-1)*2+1};
    value	= varargin{(n-1)*2+2};
    switch Parameters
        case 'Niter'
            Niter = value;
        case 'Thres'
            Thres = value;
        case 'SmoothB'
            SmoothB = value;
        case 'NStepNoLearning'
            NStepNoLearning=value;
        case 'Init'
            Init=value;
    end
end

DimX=size(X);

if iscell(Init) && length(Init)==3
    [A,B,C]=deal(Init{1},Init{2},Init{3});
else
    A=rand(DimX(1),F);
    B=rand(DimX(2),F);
    C=rand(DimX(3),F);
end

B=repmat(B,1,1,DimX(3));
Xfit=FitX(A,B,C);
Xfit_old=zeros(DimX);
iter=0;
VarX=X(:)'*X(:);

ExplaniedVariance=[];
LossOfFit=[];
Fitness=[];
Tolence=[];
%parafac parameters
ini_method=3; %initiallization with the old loading
PlotOption=0; %no plot
ScalingOption=0; %no scaling applied (hence fixed elements will not be modified
MaxIter_parafac=1;

Parafac_option=[1e-9 ini_method PlotOption ScalingOption NaN MaxIter_parafac];

parafac_const=[const(1) 0 const(3)];
I=eye(F);
parafac_OldLoad={A,I,C};
parafac_FixMode=[0 1 0];

Learn_Ratio=zeros(1,Niter);
Learn_Ratio(1:Niter-NStepNoLearning)=linspace(0,1,Niter-NStepNoLearning);
%初始点、上届，下届
while abs(Calc_Diff_Var(Xfit,Xfit_old))>Thres*VarX && iter<Niter
    iter=iter+1;
    % Step 1: Calculation of A and C
    X_reducedB=nan([DimX(1) F DimX(3)]);
    for k=1:DimX(3)
        %         PinvBk=(squeeze(B(:,:,k))*squeeze(B(:,:,k))')\squeeze(B(:,:,k));
        PinvBk=pinv(squeeze(B(:,:,k))');
        X_reducedB(:,:,k)=X(:,:,k)*PinvBk;
    end
    
    [Factors_parafac]=parafac(X_reducedB,F,Parafac_option,parafac_const,parafac_OldLoad,parafac_FixMode);
    [A,~,C]=deal(Factors_parafac{1},Factors_parafac{2},Factors_parafac{3});
    parafac_OldLoad={A,I,C};
    % Step 2:  calculation B
    %X2_31=reshape(permute(X,[2 3 1]),dimX(2),[]);
    %Bini=X2_31*pinv(krb(A,C)');
    Bini=nan(DimX(2),F);
    for j=1:DimX(2)
        Bini(j,:)=diag(pinv(A)*squeeze(X(:,j,:))*pinv(C'));
    end
    
    BBB=reshape(X,DimX(1),[])';
    BReducedA=BBB*pinv(A)';
    
    learn_ratio=Learn_Ratio(iter);
    for k=1:DimX(3)
        if const(2)==2
            bk=BReducedA((k-1)*DimX(2)+1:k*DimX(2),:)*diag(1./(C(k,:)+eps));
            bk(bk<0) = (1-learn_ratio)*bk(bk<0);
            if (rem(iter,5)==1)&&(iter<Niter-5)
                B(:,:,k)=(1-learn_ratio)*Bini+learn_ratio*bk;
            else
                B(:,:,k)=(1-learn_ratio)*squeeze(B(:,:,k))+learn_ratio*bk;
            end
            
            if SmoothB
                for j=1:F
                    if snr(squeeze(B(:,j,k)))<0
                        B(:,j,k)=sgolayfilt(squeeze(B(:,j,k)),1,11);
                    end
                end
            end
        else
            Z=(diag(C(k,:)));
            ZtZ=Z'*Z;
            ZtX=Z'*(BReducedA((k-1)*DimX(2)+1:k*DimX(2),:)');
            bk=pfls(ZtZ,ZtX,DimX(2),const(2),squeeze(B(:,:,k)),false);
            bk(bk<0) = (1-learn_ratio)*bk(bk<0);
            if (rem(iter,5)==1)&&(iter<Niter-5)
                B(:,:,k)=(1-learn_ratio)*Bini+learn_ratio*bk;
            else
                B(:,:,k)=(1-learn_ratio)*squeeze(B(:,:,k))+learn_ratio*bk;
            end
        end
        
        % Regularation of B model
        for i=1:F
            ScallBki=max(squeeze(B(:,i,k)));
            B(:,i,k)=squeeze(B(:,i,k))/(ScallBki+eps);
            C(k,i)=C(k,i)*(ScallBki+eps);
        end
    end
    Tolence=[Tolence abs(Calc_Diff_Var(Xfit,Xfit_old))];
    Xfit_old=Xfit;
    Xfit=FitX(A,B,C);
    LossOfFit=[LossOfFit Calc_Diff_Var(Xfit,X)];
    ExplaniedVariance=[ExplaniedVariance VarX-LossOfFit(end)];
    Fitness=[Fitness ExplaniedVariance(end)/VarX];
end

[A,B,C,sort_var] = SortingFacotor(A,B,C);
% calculate the entropy of B model
Xentropy=0;
for k=1:DimX(3)
    for j=1:F
        Xentropy=Xentropy+entropy(squeeze(B(:,j,k)));
    end
end
%-------------------Result------------------------
Factor={A,B,C};
Results.Eclisp=toc;
Results.F=F;
Results.Factor=Factor;
Results.X=X;
Results.Xrec=Xfit;
Results.SortVar=sort_var./VarX;
Results.iter=iter;
Results.LossOfFit=LossOfFit;
Results.ExplaniedVariance=ExplaniedVariance;
Results.Fitness=Fitness;
Results.Tolence=Tolence;
Results.Xentropy=Xentropy;
end

function Var_diff=Calc_Diff_Var(Xfit,X)
VarX=X(:)'*X(:);
VarXfit=Xfit(:)'*Xfit(:);
Var_diff=VarXfit-VarX;
end

function Xfit=FitX(A,B,C)
DimA=size(A);
DimB=size(B);
DimC=size(C);
LoadingBC=zeros(DimB(1)*DimC(1),DimB(2));
for k=1:DimC(1)
    LoadingBC((k-1)*DimB(1)+1:k*DimB(1),:)=squeeze(B(:,:,k))*diag(C(k,:));
end
Xfit=A*LoadingBC';
Xfit=reshape(Xfit,DimA(1),DimB(1),DimC(1));
end

function [A,B,C,sort_var] = SortingFacotor(A,B,C)
[nshift,F]=size(C);
var=zeros(F,1);

for i=1:F
    CBi=[];
    for j=1:nshift
        CBi=[CBi; C(j,i)*squeeze(B(:,i,j))];
    end
    Xi=A(:,i)*CBi';
    var(i)=trace(Xi'*Xi);
end
[sort_var,sort_idx]=sort(var,'descend');
A=A(:,sort_idx);
C=C(:,sort_idx);
for i=1:nshift
    B(:,:,i)=squeeze(B(:,sort_idx,i));
end
end

function [coef,yfit,exitflag]=gauss_fit(x,y,ini_coef,lb,ub,niter)
options = optimoptions('particleswarm','SwarmSize',100);
options.InitialSwarmMatrix = ini_coef;
%options.UseVectorized = true;
options.Display='off';
options.MaxIterations=niter;
GaussEqn=@(a)a(1)*exp(-((x-a(2))/a(3)).^2);
CostGauss=@(a)norm(y-a(1)*exp(-((x-a(2))/a(3)).^2));
[coef,fval,exitflag] = particleswarm(CostGauss,3,lb,ub,options);
yfit=GaussEqn(coef);
end