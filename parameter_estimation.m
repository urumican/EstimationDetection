clear all
close all

% syms phi0 phi1 phi2 x
% phi0=1;
% phi1=sqrt(3)*(3*x-1);
% phi2=sqrt(5)*(6*x^2-6*x+1);
% Phi0=[phi0, phi0];
t=10;
Nrun=1000;
a0s=ones(3,Nrun);
crlb=zeros(1,10);
MSE = zeros(3, 100);
% For our estimates
arrayMSE1 = zeros(1,100);
arrayMSE2 = zeros(1,100);
arrayMSE3 = zeros(1,100);
Bias=zeros(3,100);

for N=10:10:1000
    a_ML=zeros(3,Nrun); % ML estimator vector
    snr=ones(Nrun,1);
    
    
    for r=1:Nrun
       
        x=unifrnd(0,1,[N,1]);
%         x
        n=normrnd(0,sqrt(0.1),[N,1]);
        Phi=zeros(3,3);
        PhiV=zeros(3,1);
        y_sum=0;
%         A=zeros(3,1);
        yy1 = ones(N,1);
        yy2 = sqrt(3)*(2*x-1);
%         x.^2
        yy3 = sqrt(5)*(6*(x.^2)-6*x+1);
        Phi0=[yy1, yy2, yy3];
        Phi1=Phi0'*Phi0;
        z=Phi0*a0s(:,1);
        y=z+n;
       
        Ez=z'*z;
        En=n'*n;
        snr(r)=Ez/En;
        
        
        %for j=1:N
         %   Phi0=[1; sqrt(3)*(3*x(j)-1); sqrt(5)*(6*x(j)^2-6*x(j)+1)];
          %  Phi1=Phi0*Phi0';
           % Phi=Phi+Phi1;
            %PhiV=Phi0;
            %y_sum=sum(Phi0)+n(j);
            %A=A+y_sum*PhiV;
        %end
        
        a_ML(:,r)=inv(Phi1)*Phi0'*y;
    end
    
    E_aML=(1/Nrun)*sum(a_ML,2);
    Bias(:,N/10)=E_aML-a0s(:,1);
    
    SNR(N/10)=(1/Nrun)*sum(snr);
    
    temp1 = (a_ML-a0s);
    temp2 = temp1';
    sum1 = 0;sum2 = 0;sum3 = 0;
    
    for i=1:1:Nrun
        temp2(i,:);
        temp1(:,i);
        temp2(i,:)*temp1(:,i);
        sum1 = sum1 + temp2(i,1)*temp1(1,i);
        sum2 = sum2 + temp2(i,2)*temp1(2,i);
        sum3 = sum3 + temp2(i,3)*temp1(3,i);
    end
    
    sum1 = sum1/Nrun;
    sum2 = sum2/Nrun;
    sum3 = sum3/Nrun;
    arrayMSE1(1,N/10) = sum1;
    arrayMSE2(1,N/10) = sum2;
    arrayMSE3(1,N/10) = sum3;
    %MSE(N)=mean(sum((a_ML-a0s).^2));
    %a_Mean = sum(a_ML, 2) / size(a_ML, 2);
    %MSE(:, N)= sum((a_ML- a_Mean * ones(1, size(a_ML, 2))) .^2, 2) / size(a_ML, 2);
    
    FIM=10*N*eye(3);
    %CRLB=10/(N^2);
    CRLB=inv(FIM);
    crlb(1,N/10)=trace(CRLB)/3;
    
end

axis equal
q=10:10:1000;
loglog(q,arrayMSE1,'b','LineWidth',2);hold on
loglog(q,arrayMSE2,'y','LineWidth',2);hold on
loglog(q,arrayMSE3,'g','LineWidth',2);hold on
p=10:10:1000;
loglog(p,crlb,'r','LineWidth',2);
legend('mse1','mse2','mse3','crlb-seperate');hold off

figure
axis equal
q=10:10:1000;
loglog(q,arrayMSE1+arrayMSE2+arrayMSE3,'b','LineWidth',2);hold on
p=10:10:1000;
loglog(p,crlb*3,'r','LineWidth',2);
legend('mse','crlb');hold off

figure
axis equal
q=10:10:1000;
plot(q,Bias(1,:));%hold on
%plot(q,SNR,'r');hold off
%loglog(SNR,Bias(1,:));
