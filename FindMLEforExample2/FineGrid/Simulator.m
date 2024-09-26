function Y=Simulator(x,Level,Case) 

if Case==1
    Y=SimulatorEP(x,Level);
elseif Case==2
    Y=SimulatorPDE(x,Level);
elseif Case==3
    Y=SimulatorEnv(x,Level);    
end
end
%%
function OutputForCalibration=SimulatorEP(x,Level)
% x1 -  Outdoor air supply per unit floor area (m³/s)/m2 [0.12,0.3]
% x2 -  FCU cooling water flow rate m³/s [0.0004,0.001]
% x3 - Overnight equipment heat gain (ratio of overnight equipment heat gain to peak daytime value) n/a [0,1]

XNames={'x1' 'x2' 'x3'} ;
SimMin= [ 0.12  0.0004   0];
SimMax= [  0.3   0.001   0.975];
SimGap=SimMax-SimMin;
SimInput=SimMin + SimGap.*x ; 

%%%%%%%%%%%%%%%%%%Read all the lines in the EPBasicFile.idf
AllLines = readlines("EPBasicFile.idf");
MaxLine=length(AllLines);

%%%%%%%%%%%%%%%%%%Add values of x1 to x3 to the lines that contain variables by replacing '*****'
Lines{1}=[ 3427];
Lines{2}=[ 3808 3825 3842];
Lines{3}=[ 340 365 390 425 460 495 520 545 570 960 995];
for var=1:3
    for idx=1:numel(Lines{var})
        Line=Lines{var}(idx);
        AllLines{Line}=strrep(AllLines{Line},'*****',num2str(SimInput(var),10));
    end
end

if Level==2%HF
    Number_TimePoint_EveryHour=60;  %Number_TimePoint_EveryHour  can be 1 2 3 4 5 6 10 12 15 20 30 60
elseif Level==1%LF
    Number_TimePoint_EveryHour=5;
end

AllLines{76}=strrep(AllLines{76},'*****',num2str(Number_TimePoint_EveryHour)); %!- Number of timesteps per hour for EP simulation

%%%%%%%%%%%%%%%%%%Create the EPNewFile.idf using all the lines with added values of x1 to x3.
NewFileName='EPNewFile.idf';
fileID=fopen(NewFileName,'w+');
for idxLine=1:MaxLine
    fprintf(fileID,'%s\t\n',AllLines{idxLine});
end
fclose(fileID);

%%%%%%%%%%%%%%%%%%Run EP simulations based on the EPNewFile.idf%%%%%%%%%%%%%%%%%%
% Command=['C:\Users\PearHossain\Documents\EP\energyplus      -s C   -r   -c   -p EPNewFile        -w  EPWeather.epw      EPNewFile.idf'];
Command=['C:\EnergyPlusV9-4-0\energyplus     -s C   -r   -c   -p EPNewFile        -w  EPWeather.epw      EPNewFile.idf'];
[status,cmdout]=system(Command); %Run EP simulations

%%%%%%%%%%%%%%%%%%Read standard simulator output from the output file (EPNewFile.csv)
ReadData=readmatrix(['EPNewFile.csv']);
SimTemperature=ReadData(:,[55 67 ] );   %Simulator output with different timesteps over 0:01 - 24:00 for Sensor 1 and Sensor2

%%%%%%%%%%%%%%%%%%Average temperature for each of the two-hour intervals over 0:01 - 22:00 for Sensor 1 and Sensor2
Number_TimePoint_EveryTwoHours=Number_TimePoint_EveryHour*2;

OutputForCalibration= [ ];
for Sensor=1:2
    AverageTemperature=SimTemperature(:,Sensor); %LF simulator output with different timesteps over 0:01 - 24:00 for Sensor 1 and Sensor2
    AverageTemperature=reshape(AverageTemperature,Number_TimePoint_EveryTwoHours,[]);
    AverageTemperature11Points=mean(AverageTemperature); %Average temperature for each of the two-hour intervals
    AverageTemperature11Points(end)=[]; %Remove the last one point (22:01 to 24:00) to get 11 points
    OutputForCalibration=[OutputForCalibration AverageTemperature11Points];%2*11
end
%  delete *.audit *.bnd   *.eio  *.err  *.eso  *.end  *.epJSON  *.mdd  *.mtd  *.mtr  *.mvi  *.rdd  *.audirvaudit  *.shd  *.sql  *.html  *.svg   *.rvaudit  *.rvi  *.ZsZ  *Meter*   *Table*   *SsZ*  *ZsZ* *htm
end
%%
function [Output,Temperature]=SimulatorPDE(xx,Level)
SimMin= [20   0.0143*0.8  0.2137*0.8];
SimMax= [100  0.0143*1.2  0.2137*1.2];
SimGap=SimMax-SimMin;
SimInput=SimMin + SimGap.*xx ;
A=6.212;
ConvectionCoefficient=SimInput(1);
B=SimInput(2);
D=363.56;
E=SimInput(3);
MassDensity=8170;
T0=273.15;
AmbientTemperature=1200;InitialTemperature=300;
ThermalConductivity = @(~,state) A+B*(state.u);              
SpecificHeat=@(~,state) D+E*(state.u);                                    
                

W1=0.05*0.6;
H1=0.1*0.6;

thermalmodel = createpde('thermal','transient');
r1 = [3 4 -W1 W1 W1 -W1 -H1 -H1 H1 H1];
gdm = [r1; ]';
gm = decsg(gdm,'R1',['R1']');
geometryFromEdges(thermalmodel,gm);

if Level==1
    mesh=generateMesh(thermalmodel,'Hmax',0.0125);
elseif Level==2
    mesh=generateMesh(thermalmodel,'Hmax',0.01);
end
% pdeplot(mesh)
for i=1:thermalmodel.Geometry.NumEdges
    nodes = findNodes(mesh,'region','Edge',i);
    if(all(abs(thermalmodel.Mesh.Nodes(1,nodes)+W1)<10^-6))
        LE=i;
    elseif(all(abs(thermalmodel.Mesh.Nodes(1,nodes)-W1)<10^-6))
        RE=i;
elseif(all(abs(thermalmodel.Mesh.Nodes(2,nodes)-H1)<10^-6))
        TopE=i;        
elseif(all(abs(thermalmodel.Mesh.Nodes(2,nodes)+H1)<10^-6))
        BottomE=i;
    end
end

thermalProperties(thermalmodel,'ThermalConductivity',ThermalConductivity,'MassDensity',MassDensity,'SpecificHeat',SpecificHeat);
thermalBC(thermalmodel,'Edge',BottomE,'Temperature',InitialTemperature);
if Level==1
    thermalBC(thermalmodel,'Edge',[RE LE TopE],'ConvectionCoefficient',ConvectionCoefficient,'AmbientTemperature',AmbientTemperature);
elseif Level==2
    thermalBC(thermalmodel,'Edge',[RE LE TopE],'ConvectionCoefficient',ConvectionCoefficient,'AmbientTemperature',AmbientTemperature,'Emissivity',0.05);%celsius               %@@1000+273.15;
    thermalmodel.StefanBoltzmannConstant = 5.670367e-8;                                                                                                                                                                                                 %• SI – 5.670367e-8 W/(m2·K4)
end
thermalIC(thermalmodel,InitialTemperature);

tlist=[0:90:900];
thermalresults = solve(thermalmodel,tlist);
n=11;
[xq,yq] = meshgrid(linspace(-(n-1)/n,(n-1)/n,n)*W1,linspace(-(n-1)/n,(n-1)/n,n)*H1);%11^2
Temperature=interpolateTemperature(thermalresults,xq,yq,1:numel(tlist));%121 space -10 time

Temperature(:,1)=[];
Output=Temperature(:)';

end
%%
function [Y,YRow]=SimulatorEnv(xx,Level)
% OUTPUT AND INPUTS:
% y=row vector of scaled concentrations of the pollutant at the
%     space-time vectors (s, t)
%     Its structure is:
%     y(s_1, t_1), y(s_1, t_2), ..., y(s_1, t_dt), y(s_2, t_1), ...,
%     y(s_2,t_dt), ..., y(s_ds, t_1), ..., y(s_ds, t_dt)
% xx=[M, D, L, tau]
% s=vector of locations (optional), with default value
%     [0.5, 1, 1.5, 2, 2.5]
% t=vector of times (optional), with default value
%     [0.3, 0.6, ..., 60]

% L ∈ [0.01, 3]	location of the second spill
% τ ∈ [30.01, 30.295]   	time of the second spill
% D=(0.12-0.02)*xx(1)+0.02;
% M =7+6*xx(2);
% L=(3-0.1)*xx(1)+0.1;
L=(2-0.1)*xx(1)+0.1;
tau=(34-1)*xx(2)+1;
D=0.04;
M=13;

t=[35:5:60];
s= [1:0.4:3];


ds=length(s);
dt=length(t);
dY=ds * dt;
Y=zeros(ds, dt);
% Create matrix Y, where each row corresponds to si and each column
% corresponds to tj.
for (ii=1:ds)
    si=s(ii);
    for (jj=1:dt)
        tj=t(jj);
        
        term1a=M / sqrt(4*pi*D*tj);
        if Level==1
            term1b=max(1+(-si^2 / (4*D*tj))/6,0)^6;
        else
            term1b=exp(-si^2 / (4*D*tj));
        end
        term1=term1a * term1b;
        
        term2=0;
        if (tau < tj)
            term2a=M / sqrt(4*pi*D*(tj-tau));
            if Level==1
                term2b=max(1+(-(si-L)^2 / (4*D*(tj-tau)))/6,0)^6;
            else
                term2b=exp(-(si-L)^2 / (4*D*(tj-tau)));
            end
            term2=term2a * term2b;
        end
        
        C=term1 + term2;
        Y(ii, jj)=sqrt(4*pi) * C;
    end
end

% Convert the matrix into a vector (by rows).
YRow=Y';
Y=YRow(:)';

end