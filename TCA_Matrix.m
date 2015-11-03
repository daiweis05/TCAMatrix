function TCA_Matrix()
% Written on 03/10/2008 by Wei Li at Case Western Reserve University
% Last Modified by Wei Li on Sep. 14, 2009
clc, clear
% Animal ID: RHP009
% NMR data
P.NMRdata = [0.125 0.195 0.601 0.295 0.635 0 0.309 0.690];
% GC-MS data
CITMid = [ 6.99 3.01 10.81 14.57 21.15 34.51 8.97];
aKGMid = [ 15.80 2.03 11.48 13.30 20.91 36.49];
MALMid = [ 18.42 7.18 14.82 20.26 30.56];
SUCMID = [ 27.20 7.18 14.82 20.26 30.56];
P.CIT_Mid=CITMid/sum(CITMid);
P.aKG_Mid=aKGMid/sum(aKGMid);
P.SUC_Mid=SUCMID/sum(SUCMID);
P.MAL_Mid=MALMid/sum(MALMid);
P.CIT_Rel=P.CIT_Mid(2:end)/sum(P.CIT_Mid(2:end));
P.aKG_Rel=P.aKG_Mid(2:end)/sum(P.aKG_Mid(2:end));
P.SUC_Rel=P.SUC_Mid(2:end)/sum(P.SUC_Mid(2:end));
P.MAL_Rel=P.MAL_Mid(2:end)/sum(P.MAL_Mid(2:end));

%% Method Selection
method=4;
switch method
 case 1 % PID
 P.GCMSWT=0;
 P.GCMSRelWT=0;
 P.NMRWT=1;
 case 2 % Absolute MID
 P.GCMSWT=1;
 P.GCMSRelWT=0;
 P.NMRWT=0;
 case 3 % Relative MID
 P.GCMSWT=0;
 P.GCMSRelWT=1;
 P.NMRWT=0;
 case 4 % PID + Absolute MID
 P.GCMSWT=1;
 P.GCMSRelWT=0;
 P.NMRWT=1;
 case 5 % PID + Relative MID
 P.GCMSWT=0;
 P.GCMSRelWT=1;
 P.NMRWT=1;
end 

%% optimization process
global OAA
OAA=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
AnaOAA=zeros(4,4);
%AnaOAA(12)=1; %PYR123=1;
AnaOAA(1)=1; %PYR0=1;
P.AnaOAA=AnaOAA;
P.InactativeSUC=0;
yps=0;
ypc=0;
FAT120=0.8;
x0=[FAT120 yps ypc];
[x,Residual]=fminsearch(@NMRParameter,x0,optimset('TolX',1e-9),P);
[E,NMRdata,allPID,MidInfo,OAAarray]=NMRParameter(x,P);
P.FAT12=x(1);
P.Anaplerosis=x(2)+x(3);
disp(P)
%% plotting
NMRdata=[P.NMRdata' NMRdata];
figure(1)
subplot(2,2,1)
bar(NMRdata,1)
set(gca,'XTickLabel',{'C2D23';'C2D12';'C2Q';'C3D';'C3T';'C4D34'; ...
'C4D45'; 'C4Q'})
legend('Experiment','Simulation')
title('NMR: RHP009')
ylim([0,1])
subplot(2,3,3)
bar([P.aKG_Mid; MidInfo.aKGMID]',1)
xlim([0,7])
title('GC-MS: a-ketoglutarate')
set(gca,'XTickLabel',{'M0';'M1';'M2';'M3';'M4';'M5'})
subplot(2,3,4)
bar([P.SUC_Mid; MidInfo.SUCMID]',1)
title('GC-MS: Succinate')
xlim([0,6])
set(gca,'XTickLabel',{'M0';'M1';'M2';'M3';'M4'})
subplot(2,3,5)
bar([P.MAL_Mid; MidInfo.MALMID]',1)
xlim([0,6])
title('GC-MS: Malate')
set(gca,'XTickLabel',{'M0';'M1';'M2';'M3';'M4'})
subplot(2,3,6)
bar([P.CIT_Mid; MidInfo.CITMID]',1) 

title('GC-MS: Citrate')
xlim([0,8])
set(gca,'XTickLabel',{'M0';'M1';'M2';'M3';'M4';'M5';'M6'})
colormap gray
%% main subroutine
%----------------------------------------------------------------------
function [E,NMRdata,AllPID,MidInfo,OAAarray]=NMRParameter(x0,P)
FAT12=x0(1);
yps=x0(2);
ypc=x0(3);
InactativeSUC=P.InactativeSUC;
NaturalAbundance=0.011;
OneCarbon=[1-NaturalAbundance, NaturalAbundance];
TwoCarbon=OneCarbon'*OneCarbon;
TwoCarbon=TwoCarbon(:)';
AnaSUC=TwoCarbon'*TwoCarbon;
AcCoA=FAT12*[0 0 0 1]+(1-FAT12)*TwoCarbon;
AnaSUC_GCMS=zeros(4,4);
AnaSUC_GCMS(1,1)=1;
AcCoA_GCMS=FAT12*[0 0 0 1]+(1-FAT12)*[1 0 0 0];
%% Iterative operation
global OAA
OAA=OAA/sum(OAA);
i=0;
while i<500;
 UpdatedOAA=TCACycle(OAA,AcCoA,yps,ypc,AnaSUC,P.AnaOAA);
 i=i+1;
 OAAarray(1,:)=OAA;
 OAA=UpdatedOAA;
 OAAarray(2,:)=OAA;
 Difference=sum(abs(OAAarray(2,:)-OAAarray(1,:)));
 if Difference<1e-7;
 break
 end
end
[Temp0,Temp1,GLUT]=TCACycle(OAA,AcCoA,yps,ypc,AnaSUC,P.AnaOAA);
i=0;
while i<500;
 UpdatedOAA=TCACycle(OAA,AcCoA_GCMS,yps,ypc,AnaSUC_GCMS,P.AnaOAA);
 i=i+1;
 OAAarray(1,:)=OAA;
 OAA=UpdatedOAA;
 OAAarray(2,:)=OAA;
 Difference=sum(abs(OAAarray(2,:)-OAAarray(1,:)));
 if Difference<1e-7;
 break
 end
end
[OAA,CIT,aKG,SUC]=TCACycle(OAA,AcCoA_GCMS,yps,ypc,AnaSUC_GCMS,P.AnaOAA);
FuMal=OAA;
SUC=(SUC+AnaSUC*InactativeSUC)/(1+InactativeSUC); 

%% Calculation of MID
CITMassTemplate = [ 0 1 1 2 1 2 2 3 1 2 2 3 2 3 3 4
 1 2 2 3 2 3 3 4 2 3 3 4 3 4 4 5
 1 2 2 3 2 3 3 4 2 3 3 4 3 4 4 5
 2 3 3 4 3 4 4 5 3 4 4 5 4 5 5 6 ]';
SUCMassTemplate = [ 0 1 1 2;1 2 2 3;1 2 2 3;2 3 3 4 ];
OAAMassTemplate = [ 0 1 1 2 1 2 2 3 1 2 2 3 2 3 3 4 ]';
FuMalMassTemplate= OAAMassTemplate;
aKGMassTemplate =CITMassTemplate(1:8,:);
for i=0:6
 CITMID(i+1)=sum(sum(CIT.*(CITMassTemplate==i)));
end
for i=0:4
 SUCMID(i+1)=sum(sum(SUC.*(SUCMassTemplate==i)));
 OAAMID(i+1)=sum(sum(OAA.*(OAAMassTemplate==i)));
 MALMID(i+1)=sum(sum(FuMal.*(FuMalMassTemplate==i)));
end
for i=0:5
 aKGMID(i+1)=sum(sum(aKG.*(aKGMassTemplate==i)));
end
MidInfo.CITMID=CITMID;
MidInfo.SUCMID=SUCMID;
MidInfo.OAAMID=OAAMID;
MidInfo.aKGMID=aKGMID;
MidInfo.MALMID=MALMID;
MidInfo.MAL_Rel=MALMID(2:end)/sum(MALMID(2:end));
MidInfo.CIT_Rel=CITMID(2:end)/sum(CITMID(2:end));
MidInfo.SUC_Rel=SUCMID(2:end)/sum(SUCMID(2:end));
MidInfo.aKG_Rel=aKGMID(2:end)/sum(aKGMID(2:end));
%% Calculation of PID
C1F = sum(sum(GLUT(5:8,:)));
C2F = sum(sum(GLUT(2:2:8,:)));
C3F = sum(sum(GLUT([3 4 7 8],:)));
C4F = sum(sum(GLUT(:,3:4)));
C5F = sum(sum(GLUT(:,[2 4])));
C1S = sum(sum(GLUT([5 7],:)))/C1F;
C1D = sum(sum(GLUT([6 8],:)))/C1F;
C2S = sum(sum(GLUT(2,:)))/C2F;
C2D23=sum(sum(GLUT(4,:)))/C2F;
C2D12=sum(sum(GLUT(6,:)))/C2F;
C2Q = sum(sum(GLUT(8,:)))/C2F;
C3S = sum(sum(GLUT([3 7],1:2)))/C3F;
C3D =(sum(sum(GLUT([4 8],1:2)))+ sum(sum(GLUT([3 7],3:4))))/C3F;
C3T = sum(sum(GLUT([4 8],3:4)))/C3F;
C4S = sum(GLUT([1 2 5 6],3))/C4F;
C4D34=sum(GLUT([3 4 7 8],3))/C4F;
C4D45=sum(GLUT([1 2 5 6],4))/C4F;
C4Q = sum(GLUT([3 4 7 8],4))/C4F;
C5S = sum(GLUT(:,2))/C5F;
C5D = sum(GLUT(:,4))/C5F;
AllPID=[C2S,C2D12,C2D23,C2Q,C3S,C3D,C3T,C4S,C4D34,C4D45,C4Q]; 

NMRdata=[C2D23;C2D12;C2Q;C3D;C3T;C4D34;C4D45;C4Q];
E_NMR=sum((P.NMRdata-NMRdata').^2);
E(1)=sum((MidInfo.CITMID-P.CIT_Mid).^2);
E(2)=sum((MidInfo.aKGMID-P.aKG_Mid).^2);
E(3)=sum((MidInfo.SUCMID-P.SUC_Mid).^2);
E(4)=sum((MidInfo.MALMID-P.MAL_Mid).^2);
E_GCMS=sum(E);
E_Rel(1)=sum((MidInfo.CIT_Rel-P.CIT_Rel).^2);
E_Rel(2)=sum((MidInfo.aKG_Rel-P.aKG_Rel).^2);
E_Rel(3)=sum((MidInfo.SUC_Rel-P.SUC_Rel).^2);
E_Rel(4)=sum((MidInfo.MAL_Rel-P.MAL_Rel).^2);
E_GCMS_Rel=sum(E_Rel);
E=E_NMR*P.NMRWT+E_GCMS*P.GCMSWT+E_GCMS_Rel*P.GCMSRelWT;
return
%% Kernel
%----------------------------------------------------------------------
function [UpdatedOAA,CIT,aKG,FinalSUC]=TCACycle(OAA,AcCoA,yps,ypc,AnaSUC,AnaOAA)
index=[1 9 3 11 5 13 7 15 2 10 4 12 6 14 8 16];
CIT=OAA*AcCoA;
aKG=[eye(8) eye(8)]*CIT;
SucCoA=[eye(4) eye(4)]*aKG;
SUC=(SucCoA+AnaSUC*yps)/(1+yps);
FinalSUC=(SUC+SUC')/2;
OAAFUMAL=(FinalSUC*(1+yps)+AnaOAA*ypc)/(1+yps+ypc);
OFM=(OAAFUMAL+OAAFUMAL')/2;
OFM=OFM(:);
UpdatedOAA=OFM(index);
return 