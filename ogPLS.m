function [U,T,subU,regressionCoef] = OGPLS(xData,yData,Meta_Path,lambda,debiasing_Coef)
%% -----------------------------------
% xData is an N*M data matrix with N samples and M independent variables(metabolites)
% yData is the dependent variable with size N*1.
% Meta_Path is a M*P matrix which indicate the relation of
% metabolite and pathway, if the m metablite belong to pathway p, then the
% corresponding element is 1, otherwise is 0.
% lambda is the penalty factor
% debiasing_Coef a 1*P vector.
% U is a  M*K matrix. K is the number of latent variable of model. The kth column of U is the kth variable weighting vector
% T is a  N*K matrix. The kth column of T is the kth score vector.
% subU is P*K a cell matrix. The p row and k column is the sub-weight vector of the pth pathway by the k component. 
% regressionCoef is a M*1 vector, which is the regression coefficient of model.
%%------------------------------------------

[numSample,numMetabolite]=size(xData);
numPath=size(Meta_Path,2);
X = xData-ones(size(xData,1),1)*mean(xData);
X = X*diag(1./std(X));
Y = yData;

for g = 1:numPath
    index{g}=find(Meta_Path(:,g)==1);
    P{g}=zeros(numMetabolite,length(index{g}));
    kCoef(g)=debiasing_Coef(g); % kCoef-number of metabolite in pathways
    for i=1:length(index{g})
        P{g}(index{g}(i),i)=1;
    end
end
maxIteration = 10000;
fFunc = zeros(maxIteration,1);
fFunc(1)=1000000;
for iComp=1:min(numSample,numMetabolite)
    for g = 1:numPath
         {iComp,g}=zeros(length(index{g}),1);
    end
end
u0 = zeros(numMetabolite,1);
subU = subU0;
for iComp=1:min(numSample,numMetabolite)
    M=X'*Y;
    for iIteration = 2:maxIteration
        randPi = randperm(numPath);
        for g = 1:numPath
            uVector = u0;
            for k=1:numPath
                uVector = uVector+P{k}*subU{iComp,k};
            end
            F= uVector-P{randPi(g)}*subU{iComp,randPi(g)};
            w_soft = 2*P{randPi(g)}'*(M-F);
            norm_w_soft = norm(w_soft);
            if norm_w_soft<lambda*kCoef(randPi(g))
                subU{iComp,randPi(g)}=subU0{iComp,randPi(g)};
            else
                subU{iComp,randPi(g)} = (1-lambda*kCoef(randPi(g))/norm_w_soft)*w_soft*0.5;
            end
        end
        fFunc1=0;
        for g=1:numPath
            fFunc1 = fFunc1+(kCoef(g))*norm(subU{iComp,g});
        end
        fFunc0 = norm(M-uVector).^2;
        fFunc(iIteration) = fFunc0+lambda*fFunc1;
        if abs(1-fFunc(iIteration-1)/fFunc(iIteration))<1e-10 %| fFunc(iIteration-1)/fFunc(iIteration)<1
            break
        end
    end
    if norm(uVector)<eps
        if iComp>1
            subU(iComp:end,:)=[];
            iComp = iComp-1;
        else
            U =[];
            T=[];
            LL=[];
            iComp = iComp-1;
        end
        regression_coef=[];
        break
    end
    
    U(:,iComp)=uVector; % U is the weighting coefficent w in PLS model
    tVector=X*uVector;
    T(:,iComp)=tVector;
    normT = tVector/(tVector'*tVector);
    pVector=X'*normT;
    PP(:,iComp)=pVector;
    X=X-tVector*pVector';
    qVector=Y'*normT;
    Q(:,iComp)=qVector;
    Y=Y-tVector*qVector';
end

for i = 1:iComp
    U_star(:,i)=U(:,i);
    for j = 1:i-1
        U_star(:,i)=(eye(length(U(:,i)))-U(:,j)*PP(:,j)')*U_star(:,i);
    end
end

for i=1:iComp
    if i==1
        regression_coef(:,1)=U_star(:,1)*Q(:,1)';
    else
        regression_coef(:,i)=regression_coef(:,i-1)+U_star(:,i)*Q(:,i)';
    end
end
regressionCoef = regression_coef(:,iComp);
end
