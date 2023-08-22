clear all;
hold on;
clc;
ITER = 5000;
K = 10;
Mv = 20:30:500;
EudB = 10;
Eu = 10^(EudB/10);
rateMRC = zeros(1,length(Mv));
boundMRC = zeros(1,length(Mv));

for ite=1:ITER
    ite;
    D = Dmatrix(K);
    ubeta = diag(D);
    for mx=1:length(Mv)
        M = Mv(mx);
        pu = Eu;
%         pu = Eu/M;
%         pu = Eu/sqrt(M);
         Pp=K*pu;
         Phi=sqrt(1/K)*dftmtx(K);
         H=sqrt(1/2)*(randn(Mv(mx),K)+1i*randn(Mv(mx),K));
         G=H*sqrt(D);
         N=sqrt(1/2)*(randn(M,K)+1i*randn(M,K));
         RxBlk=sqrt(Pp)*G*Phi + N;
         Ghat=sqrt(1/Pp)*RxBlk*Phi';
         g1hat=Ghat(:,1);
         g1=G(:,1);
         e1=g1hat-g1;


        NrMRC = pu*norm(g1)^2;
        NrBoundMRC = pu*M*ubeta(1);
        DrBoundMRC = 1/K + (ubeta(1)+1/K/pu)/ubeta(1);
        DrMRC = pu*abs(g1'/norm(g1)*e1)^2 + norm(g1hat)^2/norm(g1)^2;
        DrMRC = DrMRC + pu*norm(g1hat'/norm(g1)*G(:,2:K))^2;
        DrBoundMRC = DrBoundMRC + pu*(ubeta(1)+1/K/pu)*sum(ubeta(2:K))/ubeta(1);
        rateMRC(mx) = rateMRC(mx)+ log2(1+NrMRC/DrMRC);
        boundMRC(mx) = boundMRC(mx)+ log2(1+NrBoundMRC/DrBoundMRC);
                
        
    end
end
rateMRC = rateMRC/ITER;
boundMRC = boundMRC/ITER;

figure;
plot(Mv,rateMRC,'b - s','LineWidth',3,'MarkerFaceColor','blue','MarkerSize',8.0)
hold on
plot(Mv,boundMRC,'g -. ','LineWidth',3,'MarkerFaceColor','green','MarkerSize',8.0)
grid on
title('Information Rate of Massive MIMO with Channel Estimation')
legend('MRC', 'Bound MRC','Location','SouthEast');
xlabel('Number of BS Antennas')
ylabel('Uplink Sum Rate (bits/s/Hz)')