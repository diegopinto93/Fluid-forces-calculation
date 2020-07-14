clear;
clc;
TVDm= xlsread('BPF.xlsx', 'E4:E428');%Total vertical Depth [m]
TVDft= xlsread('BPF.xlsx', 'F4:F428');%Total vertical Depth [ft]
tgrad= xlsread('BPF.xlsx','J4:J8');%Temperature gradient
md=xlsread('BPF.xlsx', 'A4:A428');%Measured depth [m]
mdft=xlsread('BPF.xlsx', 'B4:B428');%Measured depth [ft]
Inclin=xlsread('BPF.xlsx', 'C4:C428');%Inclination 
DLS=xlsread('BPF.xlsx', 'G4:G428');%Dog Leg Severity 
z=size(md); 
temp=zeros(z);%Temperature vector
tref= 50;% Temperature of reference [C]
Ts=122;%Temperature of reference [F]
Rhopipe= 64.7;
mwref= xlsread('BPF.xlsx','N4:N4');%MW reference
Pref=14.7; %Pressure reference [psi]
IDpipe=xlsread('BPF.xlsx','V4:V10'); %Inner pipe diameter [in]
IDann=xlsread('BPF.xlsx','AE4:AE11');%Inner annuli diameter [in]
ODpipe= xlsread('BPF.xlsx','T4:T10'); %Outer pipe diameter [in]
ODp1= xlsread('BPF.xlsx','AF4:AF11');
Ai=xlsread('BPF.xlsx','AH4:AH10');%[m2]
Ao=xlsread('BPF.xlsx','AG4:AG11');%[m2]
Lpipe= xlsread('BPF.xlsx','R28:R32');
Lpipe2=xlsread('BPF.xlsx','AD4:AD10');%Pipe Length  [m]
Lpipe1=xlsread('BPF.xlsx','R4:R8');% Pipe Length 1 [m]
Aipipe=xlsread('BPF.xlsx','X4:X10'); % Inner Cross sectional area pipe[in2]
wair=xlsread('BPF.xlsx','Y4:Y8');
IDtj=xlsread('BPF.xlsx','P4:P8');
BC=xlsread('BPF.xlsx','AL4:AL428');%Base case Static
BC2=xlsread('BPF.xlsx','AO4:AO428');%Base Case Dynamic
SFR=xlsread('BPF.xlsx','AR4:AR428');%Side Force RIH
SFP=xlsread('BPF.xlsx','AS4:AS428');%Side Force POOH
SFS=xlsread('BPF.xlsx','AS4:AS428');%Side Force Rotating off Bottom
Aoann=(pi.*((IDann.*0.0254).^2)./4)-Ao; % Annuli cross sectional area [m2]
AIH=zeros(z); %Area annular
Abn=0.3069;%Bit nozzles area[in2]
Ty= 6.703;% Yield stress [Pa]
n= 0.787; % Behaviour index 
k= 0.149; % Consistency index [Pa*s^n]
m= 1/n;
R= zeros(z);
R1=zeros(z);
N=size(IDpipe);
N1=size(IDann);
n1=size(Lpipe1);
na=size(z);
D=(0.0254.*IDann)-(0.0254.*ODp1);%annular diameter [m]
D1=(0.0254.*IDpipe); %inside diameter [m]
Dtj=(0.0254.*IDtj);%inside diameter tool joint [m]
q=3000/(60000); %flow [m3/s] 3000 l/min CAMBIOOOOOS
Q=3000*0.264172; % flow [gpm]CAMBIOSSSS
v=q./(0.0254^2.*Aipipe);% velocity [m/s]
vpipe=zeros(z);%velocity inside pipe
vann=zeros(z); %velocity annular
e=4.6E-05; %Pipe roughness 
b=426;
Rho= zeros(z);%Density [ppg]
Rho1=zeros(z);
P=zeros(z);% Pressure [psi]
Rho(1)= mwref;%Density 
P(1)=Pref;% Pressure at cell 1
Recp1=zeros(z);%Re pipe
Reca1=zeros(z);%Re ann 
Retj=zeros(z);%Re tool joint
deltaP=zeros(z);% Change of pressure [psi]
P2=zeros(z);
P3=zeros(z);
P4=zeros(z);
fp1=zeros(z);%Darcy friction factor pipe 
DpDlp1=zeros(z);%dpdl inside pipe [Pa/m]
Dpp1=zeros(z);%pressure loss inside pipe [Pa]
fa1=zeros(z);%Darcy friction factor annular
DpDlann=zeros(z);%dpdl annular [Pa/m]
Dpa1=zeros(z);%pressure loss inside pipe [Pa]
VFp1=zeros(z);%Viscous Force inside pipe 
VFa1=zeros(z);%Viscous Force annular
tol=0.001;% tolerance 
dfvd=zeros(z); %Additional force due to viscous drag
dptj=zeros(z);
ktj=zeros(z);
atj=zeros(z);%Area of the Tool Joint
Ao1=zeros(z);
Ai1=zeros(z);
Atj=(pi.*((IDtj).^2)./4);%CAMBIEEEE de m2 para in2
Pwem=zeros(z);%Projected weight of the mud
Pwei=zeros(z);
Pw=zeros(z);
FdeltaA= xlsread('BPF.xlsx', 'AJ4:AJ428');%HYDROSTATIC FORCES DUE TO CHANGE IN AREA
Dfpv=zeros(z); %
BHP1=zeros(z);%Effective tension Dynamic
cf=12/231; %Conversion Factor
Ptp=zeros(z); %Total pressure pipe
Pta=zeros(z); %Total pressure annular
PAF=zeros(z); %Pressure area force 
Dpa=zeros(z);
BFn=zeros(z);
BFd=zeros(z);
Ao1=([Ao1;0]);
IDH=zeros(z);%Inner diameter hole
PAFD=zeros(z);
IF=zeros(z);
IFt=zeros(z);
IFn=zeros(z);
l=zeros(z);
IFSt=zeros(z);%Inertial force tangential direction
IFSn=zeros(z);%Inertial force normal direction
IFS=zeros(z);%Inertial force total 
IF2=zeros(z);
theta=zeros(z);
PDDP=zeros(z);%Presure difference Drill Pipe
PDHWDP=zeros(z);%Presure difference HWDP
PDJar=zeros(z);%Presure difference Jar
PDDC=zeros(z);%Presure difference Drill Collar
DragR=zeros(z);%Drag force for RIH
DragP=zeros(z);%Drag force for POOH
DragS=zeros(z);%Drag force for Rotating off Bottom
%calculation Re (pipe, ann) 

for i=2:z
 if TVDm(i)<110
     temp(i)= (tref+(TVDm(i)*tgrad(1)))*(9/5)+ 32;
 else
     if TVDm(i)<1390
      temp(i)= (tref+(TVDm(i)*tgrad(2)))*(9/5)+ 32;
     else
         if TVDm(i)<1405
       temp(i)= (tref+(TVDm(i)*tgrad(3)))*(9/5)+ 32;
         else
             if TVDm(i)<1550
        temp(i)= (tref+(TVDm(i)*tgrad(4)))*(9/5)+ 32;
             else
        temp(i)= (tref+(TVDm(i)*tgrad(5)))*(9/5)+ 32;
             end
         end
     end
 end 

    P(i)=0.052*Rho(i-1)*TVDft(i);%initial guess of pressure
    deltaP(i)=P(i);
while(deltaP(i) > tol)
    P2(i)=P(i);
    Rho(i)=roliq(P(i),temp(i));%recalculate density
    P(i)=0.052*Rho(i)*TVDft(i);%recalculate pressure
    deltaP(i)=abs(P(i)-P2(i));
end

  if md(i)<= 3987.25 
       R(i)=IDpipe(5)/2;
       Ao1(i)=(pi*(ODpipe(5))^2)/4;
       Ai1(i)=(pi*(IDpipe(5))^2)/4;
       vpipe(i)=q/(0.0254^2.*Aipipe(5));
       Pwei(i)=(wair(5)*(md(i)-md(i-1))*cosd(mean([Inclin(i),Inclin(i-1)])))*2.2046226218488;
       atj(i)=(Atj(5));%Area tool joint [in2]
       Recp1(i)=((Rho(i)/0.0083)*((vpipe(i))^(2-n))*D1(5)^(n))/((Ty/8)*(D1(5)/vpipe(i))^(n) + k*((3*(n*k*(8*vpipe(i)/D1(5))^(n)/(Ty + k*(8*vpipe(i)/D1(5))^(n)))+1)/(4*(n*k*(8*vpipe(i)/D1(5))^(n)/(Ty + k*(8*vpipe(i)/D1(5))^(n)))))^(n) *8^(n-1));
       Pwei(i)=(wair(5)*(md(i)-md(i-1))*cosd(mean([Inclin(i),Inclin(i-1)])))*2.2046226218488; %Projected weight pipe 
       if Recp1(i)<1000
           ktj(i)=0;
       else
           if Recp1(i)<=3000
           ktj(i)=1.91*log10(Recp1(i))-5.64;
           else
               if Recp1(i)<=13000
                   ktj(i)=4.66-(1.05*log10(Recp1(i)));
               else
                   ktj(i)=0.33;
               end
           end
       
       end
       dptj(i)=((Rho(i)/0.0083)*ktj(i)*(vpipe(i))^2)/2;
       if Recp1(i)>=3302.45 %turbulent inside pipe
        syms x
        eqn= (4/n)*log10(Recp1(i)*x^(1-n/2)+(e/(3.7*D1(5))))-(0.4*n)-1/(x^0.5)==0;
        fp1(i)=vpasolve(eqn,x); 
       else %Laminar inside pipe
         if Recp1(i)<= 2344.95
        fp1(i)=64/Recp1(i);
        
         end
         fp1(i)=0.3164/(Recp1(i))^(1/4);
       end 
      DpDlp1(i)=(fp1(i)*(Rho(i)/0.0083)*(vpipe(i))^(2))/(2*0.0254*IDpipe(5));
      Dpp1(i)=DpDlp1(i)*(md(i)-md(i-1));
      wbp(i)=(9.8*wair(5)+((Rho(i)*119.83)*(Ai(5)-Ao(5)))*9.8)*(md(i)-md(i-1));
  else
      
    if md(i)<=4205.25
       R(i)=IDpipe(4)/2;
     vpipe(i)=q/(0.0254^2.*Aipipe(4));
       Ao1(i)=(pi*(ODpipe(4))^2)/4;
       Ai1(i)=(pi*(IDpipe(4))^2)/4;
       Pwei(i)=(wair(4)*(md(i)-md(i-1))*cosd(mean([Inclin(i),Inclin(i-1)])))*2.2046226218488;
        atj(i)=(Atj(4));
       Recp1(i)=((Rho(i)/0.0083)*((vpipe(i))^(2-n))*D1(4)^(n))/((Ty/8)*(D1(4)/vpipe(i))^(n) + k*((3*(n*k*(8*vpipe(i)/D1(4))^(n)/(Ty + k*(8*vpipe(i)/D1(4))^(n)))+1)/(4*(n*k*(8*vpipe(i)/D1(4))^(n)/(Ty + k*(8*vpipe(i)/D1(4))^(n)))))^(n) *8^(n-1));
       if Recp1(i)<1000
           ktj(i)=0;
       else
           if Recp1(i)<=3000
           ktj(i)=1.91*log10(Recp1(i))-5.64;
           else
               if Recp1(i)<=13000
                   ktj(i)=4.66-(1.05*log10(Recp1(i)));
               else
                   ktj(i)=0.33;
               end
           end
       
       end
       dptj(i)=((Rho(i)/0.0083)*ktj(i)*(vpipe(i))^2)/2;
       if Recp1(i)>=3302.45 %turbulent inside pipe
       syms x  
      eqn= (4/n)*log10(Recp1(i)*x^(1-n/2)+(e/(3.7*D1(4))))-(0.4*n)-1/(x^0.5)==0;
        fp1(i)=vpasolve(eqn,x);
         
       else
          if Recp1(i)<= 2344.95%Laminar inside pipe
        fp1(i)=64/Recp1(i);
       
          end
          fp1(i)=0.3164/(Recp1(i))^(1/4);
       end
       DpDlp1(i)=(fp1(i)*(Rho(i)/0.0083)*(vpipe(i))^(2))/(2*0.0254*IDpipe(4));
       Dpp1(i)=DpDlp1(i)*(md(i)-md(i-1)); 

    else
           if md(i)<= 4218.5
       R(i)=IDpipe(3)/2;
      vpipe(i)=q/(0.0254^2.*Aipipe(3));
        Ao1(i)=(pi*(ODpipe(3))^2)/4;
       Ai1(i)=(pi*(IDpipe(3))^2)/4;
       Pwei(i)=(wair(3)*(md(i)-md(i-1))*cosd(mean([Inclin(i),Inclin(i-1)])))*2.2046226218488 ;
       atj(i)=(Atj(3));
       Recp1(i)=((Rho(i)/0.0083)*((vpipe(i))^(2-n))*D1(3)^(n))/((Ty/8)*(D1(3)/vpipe(i))^(n) + k*((3*(n*k*(8*vpipe(i)/D1(3))^(n)/(Ty + k*(8*vpipe(i)/D1(3))^(n)))+1)/(4*(n*k*(8*vpipe(i)/D1(3))^(n)/(Ty + k*(8*vpipe(i)/D1(3))^(n)))))^(n) *8^(n-1));
        if Recp1(i)<1000
           ktj(i)=0;
       else
           if Recp1(i)<=3000
           ktj(i)=1.91*log10(Recp1(i))-5.64;
           else
               if Recp1(i)<=13000
                   ktj(i)=4.66-(1.05*log10(Recp1(i)));
               else
                   ktj(i)=0.33;
               end
           end
       
       end
       dptj(i)=((Rho(i)/0.0083)*ktj(i)*(vpipe(i))^2)/2;
           if Recp1(i)>=3302.45 %turbulent inside pipe
       syms x
       eqn= (4/n)*log10(Recp1(i)*x^(1-n/2)+(e/(3.7*D1(3))))-(0.4*n)-1/(x^0.5)==0;
       fp1(i)=vpasolve(eqn,x);
         
       else %Laminar inside pipe
           if Recp1(i)<= 2344.95
        fp1(i)=64/Recp1(i);
         
           end
           fp1(i)=0.3164/(Recp1(i))^(1/4);
           end 
        DpDlp1(i)=(fp1(i)*(Rho(i)/0.0083)*(vpipe(i))^(2))/(2*0.0254*IDpipe(3));
        Dpp1(i)=DpDlp1(i)*(md(i)-md(i-1));
      
          else
          if md(i)<= 4238.17
        Ao1(i)=(pi*(ODpipe(2))^2)/4;
       Ai1(i)=(pi*(IDpipe(2))^2)/4;
        R(i)=IDpipe(2)/2;
        vpipe(i)=q/(0.0254^2.*Aipipe(2));
        Pwei(i)=(wair(2)*(md(i)-md(i-1))*cosd(mean([Inclin(i),Inclin(i-1)])))*2.2046226218488;
        atj(i)=(Atj(2));
        Recp1(i)=((Rho(i)/0.0083)*((vpipe(i))^(2-n))*D1(2)^(n))/((Ty/8)*(D1(2)/vpipe(i))^(n) + k*((3*(n*k*(8*vpipe(i)/D1(2))^(n)/(Ty + k*(8*vpipe(i)/D1(2))^(n)))+1)/(4*(n*k*(8*vpipe(i)/D1(2))^(n)/(Ty + k*(8*vpipe(i)/D1(2))^(n)))))^(n) *8^(n-1));
          if Recp1(i)<1000
           ktj(i)=0;
       else
           if Recp1(i)<=3000
           ktj(i)=1.91*log10(Recp1(i))-5.64;
           else
               if Recp1(i)<=13000
                   ktj(i)=4.66-(1.05*log10(Recp1(i)));
               else
                   ktj(i)=0.33;
               end
           end
       
       end 
        dptj(i)=((Rho(i)/0.0083)*ktj(i)*(vpipe(i))^2)/2;
        if Recp1(i)>=3302.45 %turbulent inside pipe
       syms x
       
       eqn= (4/n)*log10(Recp1(i)*x^(1-n/2)+(e/(3.7*D1(2))))-(0.4*n)-1/(x^0.5)==0;
        fp1(i)=vpasolve(eqn,x);
       
       else %Laminar inside pipe
           if Recp1(i)<= 2344.95
        fp1(i)=64/Recp1(i);
            
           end
         fp1(i)=0.3164/(Recp1(i))^(1/4);
           end
        DpDlp1(i)=(fp1(i)*(Rho(i)/0.0083)*(vpipe(i))^(2))/(2*0.0254*IDpipe(2));
        Dpp1(i)=DpDlp1(i)*(md(i)-md(i-1));
       
          else
        R(i)=IDpipe(1)/2;
       vpipe(i)=q/(0.0254^2.*Aipipe(1));
       Ao1(i)=(pi*(ODpipe(1))^2)/4;
       Ai1(i)=(pi*(IDpipe(1))^2)/4;
       Pwei(i)=(wair(1)*(md(i)-md(i-1))*cosd(mean([Inclin(i),Inclin(i-1)])))*2.2046226218488;
         atj(i)=(Atj(1));
         Recp1(i)=((Rho(i)/0.0083)*((vpipe(i))^(2-n))*D1(1)^(n))/((Ty/8)*(D1(1)/vpipe(i))^(n) + k*((3*(n*k*(8*vpipe(i)/D1(1))^(n)/(Ty + k*(8*vpipe(i)/D1(1))^(n)))+1)/(4*(n*k*(8*vpipe(i)/D1(1))^(n)/(Ty + k*(8*vpipe(i)/D1(1))^(n)))))^(n) *8^(n-1));
     
        if Recp1(i)>=3302.45 %turbulent inside pipe
       syms x
      
       eqn= (4/n)*log10(Recp1(i)*x^(1-n/2)+(e/(3.7*D1(1))))-(0.4*n)-1/(x^0.5)==0;
        fp1(i)=vpasolve(eqn,x);
        
       else %Laminar inside pipe
           if Recp1(i)<= 2344.95
        fp1(i)=64/Recp1(i);
        
           end
           fp1(i)=0.3164/(Recp1(i))^(1/4);
        end
        DpDlp1(i)=(fp1(i)*(Rho(i)/0.0083)*(vpipe(i))^(2))/(2*0.0254*IDpipe(1));
        Dpp1(i)=DpDlp1(i)*(md(i)-md(i-1));
        
           end
            end
    end
 
  end 
   VFp1(i)=(Dpp1(i)*(pi*(R(i)*0.0254)^2))*0.2248090795; %Viscous force in the wall of the pipe [Pound]
    
   if md(i)<=180
     R1(i)=(IDann(8)-ODp1(8))/2;
     AIH(i)=((pi*(IDann(8))^2)/4)-Ao1(i);%Area inside the hole [in2]
     vann(i)=q/(Aoann(8));
     Reca1(i)=((Rho(i)/0.0083)*((vann(i))^(2-n))*D(8)^(n))/((Ty/8)*(D(8)/vann(i))^(n) + k*((3*(n*k*(8*vann(i)/D(8))^(n)/(Ty + k*(8*vann(i)/D(8))^(n)))+1)/(4*(n*k*(8*vann(i)/D(8))^(n)/(Ty + k*(8*vann(i)/D(8))^(n)))))^(n) *8^(n-1));
      if Reca1(i)>=3302.45 %turbulent annuli
       syms x
        eqn3= -2*log10((2.51/((x^0.5)*Reca1(i)))+(e/(3.7*D(8))))-1/(x^0.5)==0;
        fa1(i)=vpasolve(eqn3,x);
        
        else %Laminar annuli 
        if Reca1(i)<= 2344.95
         fa1(i)=64/Reca1(i);
         DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(8)-ODp1(8)));
         Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));
         
        end
        fa1(i)=0.3164/(Reca1(i))^(1/4); %Intermediate flow
    end
    DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(8)-ODp1(8)));
    Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));
    dfvd(i)=(pi*Dpa1(i)*((0.0254*IDann(8))^(2) - (0.0254*ODp1(8))^(2))*0.0254*ODp1(8))/(4*0.0254*(IDann(8)-ODp1(8)));
    else
    if md(i)<=1879 %Annular
     AIH(i)=((pi*(IDann(7))^2)/4)-Ao1(i);%Area inside the hole [in2]
     R1(i)=(IDann(7)-ODp1(7))/2;
     vann(i)=q/(Aoann(7));
     Reca1(i)=((Rho(i)/0.0083)*((vann(i))^(2-n))*D(7)^(n))/((Ty/8)*(D(7)/vann(i))^(n) + k*((3*(n*k*(8*vann(i)/D(7))^(n)/(Ty + k*(8*vann(i)/D(7))^(n)))+1)/(4*(n*k*(8*vann(i)/D(7))^(n)/(Ty + k*(8*vann(i)/D(7))^(n)))))^(n) *8^(n-1));
     
    if Reca1(i)>=3302.45 %turbulent annuli
       syms x
        eqn3= -2*log10((2.51/((x^0.5)*Reca1(i)))+(e/(3.7*D(7))))-1/(x^0.5)==0;
        fa1(i)=vpasolve(eqn3,x);
        
    else %Laminar annuli 
        if Reca1(i)<= 2344.95
         fa1(i)=64/Reca1(i);
         DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(7)-ODp1(7)));
         Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));
         
        end
        fa1(i)=0.3164/(Reca1(i))^(1/4); %Intermediate flow
    end
    DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(7)-ODp1(7)));
    Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));
    dfvd(i)=(pi*Dpa1(i)*((0.0254*IDann(7))^(2) - (0.0254*ODp1(7))^(2))*0.0254*ODp1(7))/(4*0.0254*(IDann(7)-ODp1(7)));
    else
        if md(i)<= 2917
            R1(i)=(IDann(6)-ODp1(6))/2;
            AIH(i)=((pi*(IDann(6))^2)/4)-Ao1(i);%Area inside the hole [in2]
            vann(i)=q/(Aoann(6));
            Reca1(i)=((Rho(i)/0.0083)*((vann(i))^(2-n))*D(6)^(n))/((Ty/8)*(D(6)/vann(i))^(n) + k*((3*(n*k*(8*vann(i)/D(6))^(n)/(Ty + k*(8*vann(i)/D(6))^(n)))+1)/(4*(n*k*(8*vann(i)/D(6))^(n)/(Ty + k*(8*vann(i)/D(6))^(n)))))^(n) *8^(n-1));
            
        if Reca1(i)>=3302.45 %turbulent annuli
       syms x
        eqn3= -2*log10((2.51/((x^0.5)*Reca1(i)))+(e/(3.7*D(6))))-1/(x^0.5)==0;
        fa1(i)=vpasolve(eqn3,x);
        
    else %Laminar annuli 
          if Reca1(i)<= 2344.95
         fa1(i)=64/Reca1(i);
         DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(6)-ODp1(6)));
         Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));       
          end
          fa1(i)=0.3164/(Reca1(i))^(1/4); %Intermediate flow
        end
        DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(6)-ODp1(6)));
        Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));
        dfvd(i)=(pi*Dpa1(i)*((0.0254*IDann(6))^(2) - (0.0254*ODp1(6))^(2))*0.0254*ODp1(6))/(4*0.0254*(IDann(6)-ODp1(6)));
        else
            if md(i)<3987.25
            R1(i)=(IDann(5)-ODp1(5))/2;
            AIH(i)=((pi*(IDann(5))^2)/4)-Ao1(i);%Area inside the hole [in2]
            vann(i)=q/(Aoann(5));
            Reca1(i)=((Rho(i)/0.0083)*((vann(i))^(2-n))*D(5)^(n))/((Ty/8)*(D(5)/vann(i))^(n) + k*((3*(n*k*(8*vann(i)/D(5))^(n)/(Ty + k*(8*vann(i)/D(5))^(n)))+1)/(4*(n*k*(8*vann(i)/D(5))^(n)/(Ty + k*(8*vann(i)/D(5))^(n)))))^(n) *8^(n-1));
           
            if Reca1(i)>=3302.45 %turbulent annuli
       syms x
        eqn3= -2*log10((2.51/((x^0.5)*Reca1(i)))+(e/(3.7*D(5))))-1/(x^0.5)==0;
        fa1(i)=vpasolve(eqn3,x);
        
            else %Laminar annuli 
                  if Reca1(i)<= 2344.95
         fa1(i)=64/Reca1(i);
          DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(5)-ODp1(5)));
           Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));
                 end
              fa1(i)=0.3164/(Reca1(i))^(1/4);
            end
           DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(5)-ODp1(5)));
           Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));
            dfvd(i)=(pi*Dpa1(i)*((0.0254*IDann(5))^(2) - (0.0254*ODp1(5))^(2))*0.0254*ODp1(5))/(4*0.0254*(IDann(5)-ODp1(5)));
            else
               if md(i)<= 4205.25
              R1(i)=(IDann(4)-ODp1(4))/2;
              AIH(i)=((pi*(IDann(4))^2)/4)-Ao1(i);%Area inside the hole [in2]
             vann(i)=q/(Aoann(4));
             Reca1(i)=((Rho(i)/0.0083)*((vann(i))^(2-n))*D(4)^(n))/((Ty/8)*(D(4)/vann(i))^(n) + k*((3*(n*k*(8*vann(i)/D(4))^(n)/(Ty + k*(8*vann(i)/D(4))^(n)))+1)/(4*(n*k*(8*vann(i)/D(4))^(n)/(Ty + k*(8*vann(i)/D(4))^(n)))))^(n) *8^(n-1));
             
                if Reca1(i)>=3302.45 %turbulent annuli
                syms x
                eqn3= -2*log10((2.51/((x^0.5)*Reca1(i)))+(e/(3.7*D(4))))-1/(x^0.5)==0;
                fa1(i)=vpasolve(eqn3,x);
                
            else %Laminar annuli 
                  if Reca1(i)<= 2344.95
                 fa1(i)=64/Reca1(i);
                 DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(4)-ODp1(4)));
                Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));
                   end
                  fa1(i)=0.3164/(Reca1(i))^(1/4);
                end  
                DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(4)-ODp1(4)));
                Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));
                dfvd(i)=(pi*Dpa1(i)*((0.0254*IDann(4))^(2) - (0.0254*ODp1(4))^(2))*0.0254*ODp1(4))/(4*0.0254*(IDann(4)-ODp1(4)));
               else
                    if md(i)<= 4218.5
                       AIH(i)=((pi*(IDann(3))^2)/4)-Ao1(i);%Area inside the hole [in2]
                        R1(i)=(IDann(3)-ODp1(3))/2;
                        vann(i)=q/(Aoann(3));
                        Reca1(i)=((Rho(i)/0.0083)*((vann(i))^(2-n))*D(3)^(n))/((Ty/8)*(D(3)/vann(i))^(n) + k*((3*(n*k*(8*vann(i)/D(3))^(n)/(Ty + k*(8*vann(i)/D(3))^(n)))+1)/(4*(n*k*(8*vann(i)/D(3))^(n)/(Ty + k*(8*vann(i)/D(3))^(n)))))^(n) *8^(n-1));
                       
                    if Reca1(i)>=3302.45 %turbulent annuli
                       syms x
                       eqn3= -2*log10((2.51/((x^0.5)*Reca1(i)))+(e/(3.7*D(3))))-1/(x^0.5)==0;
                        fa1(i)=vpasolve(eqn3,x);
                        
                        else %Laminar annuli 
                          if Reca1(i)<= 2344.95
                        fa1(i)=64/Reca1(i);
                        DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(3)-ODp1(3)));
                        Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));
                        
                         end
                        fa1(i)=0.3164/(Reca1(i))^(1/4);
                    end
                    DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(3)-ODp1(3)));
                    Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));
                 
                    dfvd(i)=(pi*Dpa1(i)*((0.0254*IDann(3))^(2) - (0.0254*ODp1(3))^(2))*0.0254*ODp1(3))/(4*0.0254*(IDann(3)-ODp1(3)));
                    else
                        if md(i)<= 4238.17
                            R1(i)=(IDann(2)-ODp1(2))/2;
                            AIH(i)=((pi*(IDann(2))^2)/4)-Ao1(i);%Area inside the hole [in2]
                            vann(i)=q/(Aoann(2));
                            Reca1(i)=((Rho(i)/0.0083)*((vann(i))^(2-n))*D(2)^(n))/((Ty/8)*(D(2)/vann(i))^(n) + k*((3*(n*k*(8*vann(i)/D(2))^(n)/(Ty + k*(8*vann(i)/D(2))^(n)))+1)/(4*(n*k*(8*vann(i)/D(2))^(n)/(Ty + k*(8*vann(i)/D(2))^(n)))))^(n) *8^(n-1));
                            
                            if Reca1(i)>=3302.45 %turbulent annuli
                             syms x
                                eqn3= -2*log10((2.51/((x^0.5)*Reca1(i)))+(e/(3.7*D(2))))-1/(x^0.5)==0;
                                fa1(i)=vpasolve(eqn3,x);
                                
                            else %Laminar annuli 
                                  if Reca1(i)<= 2344.95
                                fa1(i)=64/Reca1(i);
                                 DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(2)-ODp1(2)));
                                    Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));
                                
                                   
                                  end
                              fa1(i)=0.3164/(Reca1(i))^(1/4);
                              
                            end
                           DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(2)-ODp1(2)));
                           Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));
                           dfvd(i)=(pi*Dpa1(i)*((0.0254*IDann(2))^(2) - (0.0254*ODp1(2))^(2))*0.0254*ODp1(2))/(4*0.0254*(IDann(2)-ODp1(2)));
                        else
                            R1(i)=(IDann(1)-ODp1(1))/2;
                           AIH(i)=((pi*(IDann(1))^2)/4)-Ao1(i);%Area inside the hole [in2]
              
                            vann(i)=q/(Aoann(1));
                            Reca1(i)=((Rho(i)/0.0083)*((vann(i))^(2-n))*D(1)^(n))/((Ty/8)*(D(1)/vann(i))^(n) + k*((3*(n*k*(8*vann(i)/D(1))^(n)/(Ty + k*(8*vann(i)/D(1))^(n)))+1)/(4*(n*k*(8*vann(i)/D(1))^(n)/(Ty + k*(8*vann(i)/D(1))^(n)))))^(n) *8^(n-1));
                           
                            if Reca1(i)>=3302.45 %turbulent annuli
                             syms x
                                eqn3= -2*log10((2.51/((x^0.5)*Reca1(i)))+(e/(3.7*D(1))))-1/(x^0.5)==0;
                                fa1(i)=vpasolve(eqn3,x);
                                
                            else %Laminar annuli 
                                  if Reca1(i)<= 2344.95
                                fa1(i)=64/Reca1(i);
                                DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(1)-ODp1(1)));
                                Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));
                                  end
                                fa1(i)=0.3164/(Reca1(i))^(1/4);
                                
                            end
                           DpDlann(i)=(fa1(i)*(Rho(i)/0.0083)*(vann(i))^(2))/(2*0.0254*(IDann(1)-ODp1(1)));
                           Dpa1(i)=DpDlann(i)*(md(i)-md(i-1));
                          dfvd(i)=(pi*Dpa1(i)*((0.0254*IDann(1))^(2) - (0.0254*ODp1(1))^(2))*0.0254*ODp1(1))/(4*0.0254*(IDann(1)-ODp1(1)));%[N]
                        end
                    end
                end
            end
        end
        end
   end
   %Buoyancy calculations
   if md(i)<= 3987.25 
    BFn(i)=(Lpipe(1))*Rho(i)*(ODpipe(5)^2 - IDpipe(5)^2);
    BFd(i)=(Rhopipe*(Lpipe(1))*(ODpipe(5)^2 - IDpipe(5)^2));
   
else
if md(i)<=4205.25
     BFn(i)=(Lpipe(2))*Rho(i)*(ODpipe(4)^2 - IDpipe(4)^2);
     BFd(i)=(Rhopipe*(Lpipe(2))*(ODpipe(4)^2 - IDpipe(4)^2));
     
else
 if md(i)<= 4218.5
     BFn(i)=(Lpipe(3))*Rho(i)*(ODpipe(3)^2 - IDpipe(3)^2);
     BFd(i)=(Rhopipe*(Lpipe(3))*(ODpipe(3)^2 - IDpipe(3)^2));
 else
if md(i)<= 4238.17
    BFn(i)=(Lpipe(4))*Rho(i)*(ODpipe(2)^2 - IDpipe(2)^2);
    BFd(i)=(Rhopipe*(Lpipe(4))*(ODpipe(2)^2 - IDpipe(2)^2));
else
    BFn(i)=(Lpipe(5))*Rho(i)*(ODpipe(1)^2 - IDpipe(1)^2);
    BFd(i)=(Rhopipe*(Lpipe(5))*(ODpipe(1)^2 - IDpipe(1)^2));
end
 end
end
end 
   
   VFa1(i)=(Dpa1(i)*(pi*(R1(i)*0.0254)^2))*0.2248090795; %Viscous force in the wall of the anullar[pound]
   Ptp(i)= P(i)+((Dpp1(i)+dptj(i))*0.000145038);%Dynamic Pressure pipe
   Pta(i)= P(i)+(Dpa1(i)*0.000145038);%Dynamic Pressure annular
   
        
  Dfpv(i)=(dfvd(i)*0.22480894244319)/(AIH(i));%Extra pressure generated due to viscous drag force ANNULAR
  xdp=cumsum(Dfpv,'reverse');
end

PlossBitNoz=((8.3e-05)*Q^(2) * Rho(425))/(Abn^(2)*0.95^(2));%pressure drop through bit nozzlez psi
Plossi= ((sum(Dpp1)+sum(dptj))*0.000145038)+ PlossBitNoz% Pressure loss pipe plus bite nozzles[psi]
Plossa=sum(Dpa1)*0.000145038% Pressure loss annuli [psi]
VFann=sum(VFa1);% Total Viscous Force ann [pound]
VFann2=VFa1;
VFpipe=sum(VFp1);% Total Viscous Force pipe [ pound]
VFpipe2=VFp1;
dfvds=sum(dfvd)*0.2248090795;%
Fdrag=(VFp1-VFa1+dfvd)*0.2248090795;%Drag force [pound]
Fbtm= P(425)*Ao1(425);%Bottom pressure force [pound]
drev=fliplr(md);
for i=2:z
    if md(i)<= 3987.25 
     PDDP(i)=(Dpa1(i))*0.000145038;
  
    else
        if md(i)<=4205.25
      PDHWDP(i)=(Dpa1(i))*0.000145038; 
        else
            if md(i)<= 4218.5
            PDJar(i)=(Dpa1(i))*0.000145038;
            else
                if md(i)<= 4238.17
                  PDDC(i)=(Dpa1(i))*0.000145038;
                else
                   
                end
            end
        end
    end
         
end
DPD=sum(PDDP);
DPHW=sum(PDHWDP);
DPJ=sum(PDJar);
DPDC=sum(PDDC);
PlossBitNoz;
BF= 1-(sum(BFn)/sum(BFd)); %Buoyancy Factor
Dpp=cumsum(Dpp1);%Pressure drop iniside pipe
Dpa=cumsum(Dpa1,'reverse');%Pressure drop annuli
W=cumsum(Pwei,'reverse');%Cumulative weight of the pipe from the Bottom to the top
  for i=z:-1:2
   Pwem(i)=Ai1(i)*(mdft(i)-mdft(i-1))*cf*Rho(i)*cosd(mean([Inclin(i),Inclin(i-1)]));
   MW=cumsum(Pwem,'reverse');
   TW=W+MW;
  end
  Pltj=(dptj.*0.000145038);%Pressure losses Tool Joint psi
  Pltjcum=cumsum(Pltj);
  W2=(cumsum(Pwei)+cumsum(Pwem))./10^3;
  PP=P(425)+(Dpp(425)*0.000145038)+(Dpa(1)*0.000145038)+PlossBitNoz+xdp(1)+Pltjcum(425);%Pump pressure [psi]
  VFtj=(Pltj.*atj);%Viscous force tool joint [pound]
  BHP=PP+P(425)-(Dpp(425)*0.000145038)-PlossBitNoz-Pltjcum(425);
  
  
  BWF=BF.*(cumsum(Pwei));%Buoyant weight 
  Rho1=(Rho).*119.83;%[Kg/m3]
  
 for i=2:z
   l(i)=md(i)-md(i-1);%Lenght of the section  
   P3(i)=BHP+((Dpa(i)*0.000145038))+xdp(i)-P(i);%Pressure profile without hyd pressure [psi]
   P4(i)=BHP+((Dpa(i)*0.000145038))+xdp(i);
   if Ao1(i)==Ao1(i+1)
   else
   PAF(i)=P3(i)*(Ao1(i+1)-Ao1(i));%Hydraulic pressure due to change of area
   PAFD(i)=P(i)*(Ao1(i+1)-Ao1(i));
   end
   
 end
 P3=([P3;P3(425)]).*6894.76;%[Pa]
P4=([P4;P4(425)]).*6894.76;%from bernoulli 
 P=([P;0]).*6894.76;
 Ai1=Ai1.*0.0254^2;%[m2]
 var=(Rho1.*q^2)./(Ai1);
 
 for i=2:z
        IFSt(i)=(Ai1(i)*(P4(i)-P4(i+1)*cosd(l(i)*DLS(i))))-(var(i)*(cosd(l(i)*DLS(i))-1));
        IFSn(i)=((Ai1(i)*P4(i+1)*sind(l(i)*DLS(i)))+var(i)*sind(l(i)*DLS(i)));
 IF(i)= sqrt(((IFSt(i))^2 + (IFSn(i))^2))*0.2248090795;
 theta(i)=atand(IFSn(i)/IFSt(i));
 DragR(i)=(SFR(i)*(md(i)-md(i-1))*0.2)/10^3;%Drag force for RIH
 DragP(i)=(SFP(i)*(md(i)-md(i-1))*0.2)/10^3;%Drag force for POOH
 DragS(i)=(SFS(i)*(md(i)-md(i-1))*0.2)/10^3;%Drag force for Rotating off Bottom 
 end
VFW=BC-(VFpipe/10^3)+(VFann/10^3)-(VFtj./10^3)-(PlossBitNoz/10^3)-DragS;%Effect of viscous force
VFM=(VFann/10^3)-(VFpipe/10^3)-(VFtj./10^3)-(PlossBitNoz/10^3)-DragS;
IFC=BC-((IF)./10^3);%Inertial force
PF1=W2+cumsum(PAF./10^3)-DragS;% Pressure area force for Rotatibg Bottom 
PAF(425)=0;
Bwei=BF*(sum(Pwei));%Buoyed weight [Lb]
B=(BF.*(cumsum(Pwei))./10^3)-DragS;%Buoyancy for Rotating off Bottom 
BD=W2+cumsum(PAFD./10^3);
TF=TW(1)+sum(PAF);
Dif=W2(425)-TF; %Difference between dynamic and static pressure area
AFT=BC+(VFM)-((IF)./10^3)+(((PAF)./10^3)-DragS); %All forces together effect
AFTR=AFT-DragR;%All forces together effect RIH
AFTP=AFT+DragP;%All forces together effect POOH
BR=B-DragR;%Effect of bouyancy force RIH
BP=B+DragP;%Effect of piston force POOH
PAR=PF1-DragR;%Effect of piston force RIH
PAP=PF1+DragP;%Effect of piston force POOH
IFR=IFC-DragR;%Effect of inertal force RIH
IFP=IFC+DragP;%Effect of inertal force POOH
VFR=VFW-DragR;%Effect of viscous force RIH
VFP=VFW+DragP;%Effect of viscous force POOH
EIF=((BC-IFC)./BC).*100;%Effect of Inertial Force
plot(BC,md,'m',IFC,md,'r',BC2,md,'k')%
set(gca, 'YDir','reverse')
title('HKLD vs MD')
xlabel('Weight [Kpound]') 
ylabel('Depth [m]') 
legend('Base Case Static','Effect of All Forces','Base Case dynamic')
grid on
%% Effect of inertial force versus inclination
plot(Inclin,md,EIF,md)%
set(gca, 'YDir','reverse')
title('Inclination vs Effect of inertial force')
xlabel('Inclination [degrees]') 
ylabel('Depth [m]') 
legend('Inclination','Effecf')
grid on
%% Pressure Area
plot(PF1,md,'r',BC2,md,'k',BC,md,'b')
set(gca, 'YDir','reverse')
title('P/A force vs MD')
xlabel('Weight [kpound]') 
ylabel('Depth [m]') 
legend('P/A force w/ flow effect','Base case Dyn','Base case Stat')
grid on
%% Inertial force
plot(IFC,md,'b',BC2,md,'k',BC,md,'r')
set(gca, 'YDir','reverse')
title('Inertial force vs MD')
xlabel('Weight [kpound]') 
ylabel('Depth [m]') 
legend('Inertial force effect','Base Case Dynamic','Base Case Static')
grid on
%% Viscous Force
plot(VFW,md,BC,md,'k',BC2,md,'r')
set(gca, 'YDir','reverse')
title('Viscous force vs MD')
xlabel('Weight [kpound]') 
ylabel('Depth [m]') 
legend('Effect of Viscous','Base Case Static','Base Case Dyn')
grid on
%% Effect all forces
plot(BC,md,AFT,md,'r',BC2,md,'k')%
set(gca, 'YDir','reverse')
title('HKLD vs MD')
xlabel('Weight [Kpound]') 
ylabel('Depth [m]') 
legend('Base Case Static','Effect of All Forces','Base Case Dynamic')
grid on
%% Buoyancy vs BC
plot(B,md,'b',BC,md,'r')
set(gca, 'YDir','reverse')
title('HKLD vs MD')
xlabel('Weight [kpound]') 
ylabel('Depth [m]') 
legend('Buoyant weight','Base Case')
grid on
%% Effect all forces td
plot(BC,md,'m',AFT,md,'r',AFTR,md,'k',AFTP,md,'b')%
set(gca, 'YDir','reverse')
title('HKLD vs MD')
xlabel('Weight [Kpound]') 
ylabel('Depth [m]') 
legend('Base Case Static','Effect of All Forces','RIH','POOH')
grid on
%% buoyancy T&D
plot(BC,md,'m',B,md,'r',BR,md,'k',BP,md,'b')%
set(gca, 'YDir','reverse')
title('HKLD vs MD')
xlabel('Weight [Kpound]') 
ylabel('Depth [m]') 
legend('Base Case Static','Buoyancy','RIH','POOH')
grid on
%% P/A T&D
plot(BC,md,'m',PF1,md,'r',PAR,md,'k',PAP,md,'b')%
set(gca, 'YDir','reverse')
title('HKLD vs MD')
xlabel('Weight [Kpound]') 
ylabel('Depth [m]') 
legend('Base Case Static','Pressure Area','RIH','POOH')
grid on
%% IF T&D
plot(BC,md,'m',IFC,md,'r',IFR,md,'k',IFP,md,'b')%
set(gca, 'YDir','reverse')
title('HKLD vs MD')
xlabel('Weight [Kpound]') 
ylabel('Depth [m]') 
legend('Base Case Static','Inertial Force','RIH','POOH')
grid on
%% Viscous T&D
plot(BC,md,'m',VFW,md,'r',VFR,md,'k',VFP,md,'b')%
set(gca, 'YDir','reverse')
title('HKLD vs MD')
xlabel('Weight [Kpound]') 
ylabel('Depth [m]') 
legend('Base Case Static','Viscous','RIH','POOH')
grid on
%%
plot(BC,md,'m',AFT,md,'r',BC2,md,'k')%
set(gca, 'YDir','reverse')
title('HKLD vs MD')
xlabel('Weight [Kpound]') 
ylabel('Depth [m]') 
legend('Base Case Static','Effect of All Forces','Base Case dynamic')
grid on
%%
plot(BC,md,'m',VFW,md,'r',BC2,md,'k')%
set(gca, 'YDir','reverse')
title('HKLD vs MD')
xlabel('Weight [Kpound]') 
ylabel('Depth [m]') 
legend('Base Case Static','Effect of viscous force','Base Case dynamic')
grid on