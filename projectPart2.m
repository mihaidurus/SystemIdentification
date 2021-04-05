close all
clear all
format long
% data loading
data=load('iddata-15.mat');
uid=data.id.u;
yid=data.id.y;
uval=data.val.u;
yval=data.val.y;
Ts=data.id.Ts;
N=length(yid);
t=0:0.045:44.955;% creating time vector to be later used in plotting the results

naMax=3;% na and nb are order of dynamics
nbMax=naMax;
nkMax=2; % the delay signal
mMax=5;% Mmax->order of the system (for mMax>=6 we will see a huge overfitting effect)

% plotting the given data set
figure(1);subplot(211);plot(t,yval);xlabel('t');ylabel('yval');
title('Validation data set');
subplot(212);plot(t,uval);xlabel('t');ylabel('uval');
figure(2);subplot(211);plot(t,yid);xlabel('t');ylabel('yid');
title('Identification data set');
subplot(212);plot(t,uid);xlabel('t');ylabel('uid');

% these "best.." matrixes will be used to store the values of: na,nb,nk,
% m and the lowest MSE for each part of our design (prediction and
% simulation)
bestPred_id=zeros(mMax,4);
bestPred_val=zeros(mMax,4);
bestSim_id=zeros(mMax,4);
bestSim_val=zeros(mMax,4);

% these additional vectors will be used to store the minimum MSE for each
% degree m for later comparison of the MSE vs Degree
MSEv_pred_id=zeros(1,mMax);
MSEv_sim_id=zeros(1,mMax);
MSEv_pred_val=zeros(1,mMax);
MSEv_sim_val=zeros(1,mMax);
count=1;% used to iterate through MSEv vectors

% We iterate through each value of m (m is the order of the polynomial regressor)
for m=1:mMax
    
    % We create 4 matrixes, 2 for prediction corresponding to each data set
    % and 2 for simulation, these will be used to store the MSEs
    MSE_pred_id=zeros(naMax*nbMax*nkMax,4);
    MSE_pred_val=zeros(naMax*nbMax*nkMax,4);
    MSE_sim_id=zeros(naMax*nbMax*nkMax,4);
    MSE_sim_val=zeros(naMax*nbMax*nkMax,4);
    row=1;% used to iterate through the rows of the MSE matrixes when na,nb or nk changes
    
    % We iterate to each imposed value for na,nb,nk (na, nb are orders of the system and nk is the delay)
    for na=1:naMax
        for nb=1:nbMax
            for nk=1:nkMax
                
                % calling ARXstructure function to compute the d vectors
                % used in the regressors
                did=ARXstructure(na,nb,nk,yid,uid,N);
                dval=ARXstructure(na,nb,nk,yval,uval,N);
                
                % calling Regressors function to compute the regressors
                Phi_id=Regressors(m,did,na,nb,N);
                Phi_val=Regressors(m,dval,na,nb,N);
                
                % applying liniar regression to obtain the parameters of the model
                theta=Phi_id\yid;
                
                % using the previous parameters in order to generate the prediction model output
                Yhatpred_val=Phi_val*theta;
                Yhatpred_id=Phi_id*theta;
                
                % we calculate the MSEs and put them in the MSE matrix from
                % where we extract the row with the lowest MSE
                MSEpred_id=0;
                for i=1:length(yid)
                    MSEpred_id=MSEpred_id+1/length(yid)*((Yhatpred_id(i,1)-yid(i,1)).^2);
                end
                
                MSE_pred_id(row,1)=na;
                MSE_pred_id(row,2)=nb;
                MSE_pred_id(row,3)=nk;
                MSE_pred_id(row,4)=MSEpred_id;
                
                MSEpred_val=0;
                for i=1:length(yval)
                    MSEpred_val=MSEpred_val+1/length(yval)*((Yhatpred_val(i,1)-yval(i,1)).^2);
                end
                
                MSE_pred_val(row,1)=na;
                MSE_pred_val(row,2)=nb;
                MSE_pred_val(row,3)=nk;
                MSE_pred_val(row,4)=MSEpred_val;
                
                % We apply the same algorithm as above for the simulation part
                PhiSim_val=RegressorsSIM(uval,m,theta,na,nb,nk,N);
                PhiSim_id=RegressorsSIM(uid,m,theta,na,nb,nk,N);
                
                Yhatsim_id=PhiSim_id*theta;
                Yhatsim_val=PhiSim_val*theta;
                
                MSEsim_val=0;
                for i=1:length(yval)
                    MSEsim_val=MSEsim_val+1/length(yval)*((Yhatsim_val(i,1)-yval(i,1)).^2);
                end
                
                MSEsim_id=0;
                for i=1:length(yid)
                    MSEsim_id=MSEsim_id+1/length(yid)*((Yhatsim_id(i,1)-yid(i,1)).^2);
                end
                
                MSE_sim_id(row,1)=na;
                MSE_sim_id(row,2)=nb;
                MSE_sim_id(row,3)=nk;
                MSE_sim_id(row,4)=MSEsim_id;
                
                MSE_sim_val(row,1)=na;
                MSE_sim_val(row,2)=nb;
                MSE_sim_val(row,3)=nk;
                MSE_sim_val(row,4)=MSEsim_val;
                row=row+1;
                
            end
        end
    end
    % delete the semicolon to see the MSEpred/sim matrixes
    MSE_pred_id;
    MSE_pred_val;
    MSE_sim_id;
    MSE_sim_val;
    
    % We extract the minimum from the last column and store the number of
    % the line which represents the order of the polynomial
    min1=MSE_pred_id(1,4);% used to store the 1st minimum
    min2=MSE_pred_val(1,4);
    min3=MSE_sim_id(1,4);
    min4=MSE_sim_val(1,4);
    
    poz1=1;% used to store the lines with the lowest MSE
    poz2=1;
    poz3=1;
    poz4=1;
    
    for i=2:naMax*nbMax*nkMax
        
        if(MSE_pred_id(i,4)<min1)
            min1=MSE_pred_id(i,4);
            poz1=i;
        end
        
        if(MSE_pred_val(i,4)<min2)
            min2=MSE_pred_val(i,4);
            poz2=i;
        end
        
        if(MSE_sim_id(i,4)<min3)
            min3=MSE_sim_id(i,4);
            poz3=i;
        end
        
        if(MSE_sim_val(i,4)<min4)
            min4=MSE_sim_val(i,4);
            poz4=i;
        end
        
    end
    % here m will iterate to the next value, so we have to save the
    % correspoding parameters for each line which contains a minimum MSE
    bestPred_id(m,1)=MSE_pred_id(poz1,1);
    bestPred_id(m,2)=MSE_pred_id(poz1,2);
    bestPred_id(m,3)=MSE_pred_id(poz1,3);
    bestPred_id(m,4)=MSE_pred_id(poz1,4);
    MSEv_pred_id(count)=MSE_pred_id(poz1,4);
    
    bestPred_val(m,1)=MSE_pred_val(poz2,1);
    bestPred_val(m,2)=MSE_pred_val(poz2,2);
    bestPred_val(m,3)=MSE_pred_val(poz2,3);
    bestPred_val(m,4)=MSE_pred_val(poz2,4);
    MSEv_pred_val(count)=MSE_pred_val(poz2,4);
    
    bestSim_id(m,1)=MSE_sim_id(poz3,1);
    bestSim_id(m,2)=MSE_sim_id(poz3,2);
    bestSim_id(m,3)=MSE_sim_id(poz3,3);
    bestSim_id(m,4)=MSE_sim_id(poz3,4);
    MSEv_sim_id(count)=MSE_sim_id(poz3,4);
    
    bestSim_val(m,1)=MSE_sim_val(poz4,1);
    bestSim_val(m,2)=MSE_sim_val(poz4,2);
    bestSim_val(m,3)=MSE_sim_val(poz4,3);
    bestSim_val(m,4)=MSE_sim_val(poz4,4);
    MSEv_sim_val(count)=MSE_sim_val(poz4,4);
    
    count=count+1;
end
% delete the semicolon to check the bestPred/Sim matrixes
bestPred_id;
bestPred_val;
bestSim_id;
bestSim_val;

% After we reach mMax rows we again find the line with the minimum MSE and
% compute the best (of the best degree) parameters from it
min1=bestPred_id(1,4);
min2=bestPred_val(1,4);
min3=bestSim_id(1,4);
min4=bestSim_val(1,4);

poz1=1;
poz2=1;
poz3=1;
poz4=1;

for i=2:mMax
    if(bestPred_id(i,4)<min1)
        min1=bestPred_id(i,4);
        poz1=i;
    end
    if(bestPred_val(i,4)<min2)
        min2=bestPred_val(i,4);
        poz2=i;
    end
    if(bestSim_id(i,4)<min3)
        min3=bestSim_id(i,4);
        poz3=i;
    end
    if(bestSim_val(i,4)<min4)
        min4=bestSim_val(i,4);
        poz4=i;
    end
end

% here we extract the values from our best matrixes and assign them to some
% variables to be later used in plotting the best fits
% even tho we computed the best paramters for identification too, we will
% use in our final comparisons only the best parameters on validation
bestna_pred_id=bestPred_id(poz1,1);
bestnb_pred_id=bestPred_id(poz1,2);
bestnk_pred_id=bestPred_id(poz1,3);
bestm_pred_id=poz1;

bestna_pred_val=bestPred_val(poz2,1);
bestnb_pred_val=bestPred_val(poz2,2);
bestnk_pred_val=bestPred_val(poz2,3);
bestm_pred_val=poz2;

bestna_sim_id=bestSim_id(poz3,1);
bestnb_sim_id=bestSim_id(poz3,2);
bestnk_sim_id=bestSim_id(poz3,3);
bestm_sim_id=poz3;

bestna_sim_val=bestSim_val(poz4,1);
bestnb_sim_val=bestSim_val(poz4,2);
bestnk_sim_val=bestSim_val(poz4,3);
bestm_sim_val=poz4;

% PREDICTION PART

% using the best parameters computed above,
% parameters which lead to the best fit on the given data
did=ARXstructure(bestna_pred_val,bestnb_pred_val,bestnk_pred_val,yid,uid,N);% calling ARXstructure function based on the best parameters in order to generate the d vector which will be used for the regressors
dval=ARXstructure(bestna_pred_val,bestnb_pred_val,bestnk_pred_val,yval,uval,N);

Phi_pred_id=Regressors(bestm_pred_val,did,bestna_pred_val,bestnb_pred_val,N);% calling Regressors to generate a regressor based on the d vector and the best parameters computed earlier
Phi_pred_val=Regressors(bestm_pred_val,dval,bestna_pred_val,bestnb_pred_val,N);

theta=Phi_pred_id\yid;% applying liniar regression to obtain the parameters of the model

Ypred_id=Phi_pred_id*theta;% using the previous parameters in order to generate the prediction model output
Ypred_val=Phi_pred_val*theta;

% plotting the resulted ouputs of our model and comparing it with the given
% identification/validation data set
figure(3);subplot(221);plot(t,yid,t,Ypred_id);% plotting the prediction output on identification
xlabel('t');
legend('yid','Ypredid');
title(['Prediction on identification data, m=' num2str(bestm_pred_val) ', na=' num2str(bestna_pred_val) ', nb=' num2str(bestnb_pred_val) ', nk=' num2str(bestnk_pred_val)]);
subplot(223);plot(t,yval,t,Ypred_val);% plotting the prediction output on validation
xlabel('t');
legend('yval','Ypredval');
title(['Prediction on validation data, m=' num2str(bestm_pred_val) ', na=' num2str(bestna_pred_val) ', nb=' num2str(bestnb_pred_val) ', nk=' num2str(bestnk_pred_val)]);

% SIMULATION PART

% we apply the same algorithm in order to find the simulation model, but
% this time without calling the ARXstructure function because it is already
% done in the function RegressorsSIM (for the simulation part the d vector is
% updated with the previous computed ouputs)
% we have to recalculate the d vector and Phi matrix to have the dimensions
% coresponding to the best simulation error and the one from the prediction
% from computing different values of na,nb,nk and m, it can be seen that
% the simulation gets better as m increases up to 5, after that value we
% see overfitting

did=ARXstructure(bestna_sim_val,bestnb_sim_val,bestnk_sim_val,yid,uid,N);% calling ARXstructure function based on best parameters in order to generate the d vector which will be used for the regressors
Phi_pred_id=Regressors(bestm_sim_val,did,bestna_sim_val,bestnb_sim_val,N);% calling Regressors to generate a regressor based on the d vector and the best parameters computed earlier
theta=Phi_pred_id\yid;% applying liniar regression to obtain the parameters of the model

Phi_sim_id=RegressorsSIM(uid,bestm_sim_val,theta,bestna_sim_val,bestnb_sim_val,bestnk_sim_val,N);
Yhat_sim_id=Phi_sim_id*theta;
Phi_sim_val=RegressorsSIM(uval,bestm_sim_val,theta,bestna_sim_val,bestnb_sim_val,bestnk_sim_val,N);
Yhat_sim_val=Phi_sim_val*theta;

% plotting the resulted ouputs of our model and comparing it with the given
% identification/validation data set
subplot(222);plot(t,Yhat_sim_id,t,yid);
xlabel('t');
legend('ysimdid','yid');
title(['Simulation on identification data, m=' num2str(bestm_sim_id)  ', na=' num2str(bestna_sim_val) ', nb=' num2str(bestnb_sim_val) ', nk=' num2str(bestnk_sim_val)]);
subplot(224);plot(t,Yhat_sim_val,t,yval);
xlabel('t');
legend('ysimval','yval');
title(['Simulation on validation data, m=' num2str(bestm_sim_val) ', na=' num2str(bestna_sim_val) ', nb=' num2str(bestnb_sim_val) ', nk=' num2str(bestnk_sim_val)]);

% plotting the lowest MSE values with their respective degree for
% prediction/simulation
figure(5);subplot(211);plot(1:mMax,MSEv_pred_id);
title(['MSEpredID= ' num2str(min(MSEv_pred_id))]);
legend('MSEpredID');ylabel('MSEpred');xlabel('DEGREE (m)')
subplot(212);plot(1:mMax,MSEv_pred_val);
title(['MSEpredVAL= ' num2str(min(MSEv_pred_val))]);
legend('MSEpredVAL');ylabel('MSEpred');xlabel('DEGREE (m)')

figure(6);subplot(211);plot(1:mMax,MSEv_sim_id);
title(['MSEsimID= ' num2str(min(MSEv_sim_id))]);
legend('MSEsimID');ylabel('MSEsim');xlabel('DEGREE (m)')
subplot(212);plot(1:mMax,MSEv_sim_val);
title(['MSEsimVAL= ' num2str(min(MSEv_sim_val))]);
legend('MSEsimVAL');ylabel('MSEsim');xlabel('DEGREE (m)')

% we recompute all the above regressors and d vectors for simulation and
% prediction to show off how the results are influenced by the values of
% na, nb, nk and m 
did=ARXstructure(1,3,1,yid,uid,N);
dval=ARXstructure(1,3,1,yid,uid,N);
Phi_pred_id=Regressors(2,did,1,3,N);
Phi_pred_val=Regressors(2,dval,1,3,N);
theta=Phi_pred_id\yid;
Ypred_id=Phi_pred_id*theta;
Ypred_val=Phi_pred_val*theta;

Phi_sim_id=RegressorsSIM(uid,2,theta,1,3,1,N);
Phi_sim_val=RegressorsSIM(uval,2,theta,1,3,1,N);
Yhat_sim_id=Phi_sim_id*theta;
Yhat_sim_val=Phi_sim_val*theta;

% plots for our model outputs using parameters different than the best ones
figure(7);
subplot(221);plot(t,yid,t,Ypred_id);% plotting the prediction output on identification
xlabel('t');
legend('yid','Ypredid');
title(['Prediction on identification data, m=' num2str(2) ', na=' num2str(1) ', nb=' num2str(3) ', nk=' num2str(1)]);
subplot(223);plot(t,yval,t,Ypred_val);% plotting the prediction output on validation
xlabel('t');
legend('yval','Ypredval');
title(['Prediction on validation data, m=' num2str(2) ', na=' num2str(1) ', nb=' num2str(3) ', nk=' num2str(1)]);

subplot(222);plot(t,Yhat_sim_id,t,yid);
xlabel('t');
legend('ysimdid','yid');
title(['Simulation on identification data, m=' num2str(2) ', na=' num2str(1) ', nb=' num2str(3) ', nk=' num2str(1)]);
subplot(224);plot(t,Yhat_sim_val,t,yval);
xlabel('t');
legend('ysimval','yval');
title(['Simulation on validation data, m=' num2str(2) ', na=' num2str(1) ', nb=' num2str(3) ', nk=' num2str(1)]);



function d = ARXstructure(na,nb,nk,y,u,N)
% generate the output part of the ARX structure
k=1;
while k<=N % we iterate through each value of the data set
    for col=1:na % for the outpart we iterate just up to na
        if (k-col > 0)% negative or zero indexes are considered to be 0
            d(k,col)=y(k-col);
        else
            d(k,col)=0;
        end
    end
    k=k+1;
end

% Generate the input part of the ARX structure
k=1;
while k<=N % we iterate through each value of the data set
    for col=na+1:na+nb % for the input part we iterate from na+1 up to na+nb
        if (k-col+na-nk+1> 0)% negative or zero indexes are considered to be 0
            d(k,col)=u(k-col+na-nk+1);% for the input u we have to take into
            % consideration the delay nk too
        else
            d(k,col)=0;
        end
    end
    k=k+1;
end
end

function PHI = Regressors(mb,d,na,nb,N)
mb=mb+1;
PHI=ones(1,1);
x=zeros(1,na+nb);
for k=1:N
    % -> we get the variables from our arx structure for generating the regressor
    % -> we will consider an array (pow()) in which we will generate all the combinations of powers for our regressor
    % -> pow() will be obtained by generating all the numbers in numerical base m+1 where m = the maximum power
    x(1,:)=d(k,:);
    % ->we contruct an array used to store the powers of thre regressor
    for a=1:na+nb
        pow(a)=0;
    end
    % -> we construct an array which will be used to store the carries when we add the numbers in m+1 base
    for a=1:na+nb+1
        r(a)=0;
    end
    r(na+nb+1)=1; % -> we start from 00..1
    count=1;
    while(r(1)==0) % -> we add +1 on numbers in m+1 base until we reach the maximu number possible on na+nb digits
        i=1; %-> i is used to iterate throught the digits of our numbers in m+1 base
        while(i<=na+nb)
            pow(i)=pow(i)+r(i+1); % -> we always add the digits of our computed numbers with their corresponding carry
            
            if ((pow(i)==-1) && (i<na+nb))
                pow(i)=0;
            end
            r(i+1)=0; % -> reset the carry which was just added
            r(na+nb+1)=1; % -> last carry will always be one as we generate the numbers in ascending order one by one 0001 0002
            
            % -> if after adding 1 a certain digit exceeds the maximum value
            % -> we need to set the next carry to one and move the index i to beginning of the number
            if (pow(i) > (mb-1))
                pow(i)=-1; % not to miss the numbers of form 1000 2000
                r(i)=1;
                i=0;
            end
            i=i+1;
        end
        s=0; % -> check the sum of powers for validation
        for l=1:na+nb
            s=s+pow(l);
        end
        % -> If sum is <= maximum power we compute the value for the
        % -> regressor using the x's and their power combination just obtained
        if (s<mb)
            value=1;
            for j=1:na+nb
                value=value*x(j)^pow(j);
            end
            PHI(k,count)=value;
            count=count+1;
        end
    end
end
end


function PHI = RegressorsSIM(u,mb,theta,na,nb,nk,N)
%Generate the output part of the ARX structure
dsim=zeros(1,na+nb);
mb=mb+1;
% Generate the input part of the ARX structure
k=1;
while k<=N
    for col=na+1:na+nb
        if (k-col+na-nk+1> 0)
            dsim(k,col)=u(k-col+na-nk+1);
        else
            dsim(k,col)=0;
        end
    end
    k=k+1;
end

PHI=ones(1,1);
x=zeros(1,na+nb);
poz=1;
for k=1:N
    % -> we get the variables from our arx structure for generating the regressor
    % -> we will consider an array (pow()) in which we will generate all the combinations of powers for our regressor
    % -> pow() will be obtained by generating all the numbers in numerical base m+1 where m = the maximum power
    x(1,:)=dsim(k,:);
    % -> we contruct an array used to store the powers of thre regressor
    for a=1:na+nb
        pow(a)=0;
    end
    % -> we construct an array which will be used to store the carries when we add the numbers in m+1 base
    for a=1:nb+na+1
        r(a)=0;
    end
    r(na+nb+1)=1; % -> we start from 00..1
    count=1;
    while(r(1)==0) % -> we add +1 on numbers in m+1 base until we reach the maximu number possible on na+nb digits
        i=1; %-> i is used to iterate throught the digits of our numbers in m+1 base
        while(i<=na+nb)
            pow(i)=pow(i)+r(i+1); % -> we always add the digits of our computed numbers with their corresponding carry
            if ((pow(i)==-1) && (i<nb+na))
                pow(i)=0;
            end
            r(i+1)=0; % -> reset the carry which was just added
            r(nb+na+1)=1; % -> last carry will always be one as we generate the numbers in ascending order one by one 0001 0002
            % -> if after adding 1 a certain digit exceeds the maximum value
            % -> we need to set the next carry to one and move the index i to beginning of the number
            if (pow(i) > (mb-1))
                pow(i)=-1; % not to miss the numbers of form 1000 2000
                r(i)=1;
                i=0;
            end
            i=i+1;
        end
        s=0; % -> check the sum of powers for validation
        for l=1:na+nb
            s=s+pow(l);
        end
        % -> If sum is <= maximum power we compute the value for the
        % -> regressor using the x's and their power combination just obtained
        if (s<mb)
            value=1;
            for j=1:na+nb
                value=value*x(j)^pow(j);
            end
            PHI(k,count)=value;
            count=count+1;
        end
    end
    % -> For simulation we generate the next set of x's based on what we computed so far (step by step)
    if(k>1)
        i=na;
        while(i~=1)
            dsim(k+1,i)=dsim(k,i-1);
            i=i-1;
        end
    end
    dsim(k+1,1)=PHI(k,:)*theta;
end
end