function y=Simulator(x,Level,Case) 
%x is the vector of calibration parameters. If Case=i, the bi-fidelity simulator for Example i is called; Level=2 gives the HF output while Level=1 gives the LF output.
if Case==1
    y=SimulatorEP(x,Level);
elseif Case==2
    y=SimulatorPDE(x,Level);
elseif Case==3
    y=SimulatorEnv(x,Level);    
end
end
%%
function OutputForCalibration=SimulatorEP(x,Level)
% Inputs(1) =  Outdoor air supply per unit floor area; unit=(m^3/s)/m^2; range=[0.12,0.3].
% Inputs(2) =  FCU cooling water flow rate; unit=m^3/s; range=[0.0004,0.001].
% Inputs(3) =  Overnight equipment heat gain (expressed as a fraction of the peak daytime equipment heat gain); unit=n/a; range=[0,0.975].

MinInputs=[0.12 0.0004 0];
MaxInputs=[0.3  0.001  0.975];
RangeInputs=MaxInputs-MinInputs;
Inputs=MinInputs+RangeInputs.*x; 

%Read all the lines in the EPBasicFile.idf file.
AllLines=readlines("EPBasicFile.idf");
MaxLine=length(AllLines);

%Replaces '*****' with the value of Inputs(i) in the lines that need the value of this input for each i = 1, 2, 3.
Lines{1}=[3427];
Lines{2}=[3808 3825 3842];
Lines{3}=[340 365 390 425 460 495 520 545 570 960 995];
for inputno=1:3
    for idx=1:numel(Lines{inputno})
        Line=Lines{inputno}(idx);
        AllLines{Line}=strrep(AllLines{Line},'*****',num2str(Inputs(inputno),10));
    end
end

if Level==2%HF
    Number_TimePoints_EveryHour=60; %Number_TimePoints_EveryHour can be 1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, or 60; it is the number of timesteps per hour for the EP simulation.
elseif Level==1%LF
    Number_TimePoints_EveryHour=5;
end

AllLines{76}=strrep(AllLines{76},'*****',num2str(Number_TimePoints_EveryHour)); %Sets the number of timesteps per hour for the EP simulation.

%Creates the EPNewFile.idf file using the AllLines variable.
NewFileName='EPNewFile.idf';
fileID=fopen(NewFileName,'w+');
for idxLine=1:MaxLine
    fprintf(fileID,'%s\t\n',AllLines{idxLine});
end
fclose(fileID);

%%%%%%%%%%%%%%%%%%Runs the EP simulation based on the EPNewFile.idf file.%%%%%%%%%%%%%%%%%%
Command=['C:\Users\PearHossain\Documents\EP\energyplus      -s C   -r   -c   -p EPNewFile        -w  EPWeather.epw      EPNewFile.idf'];
%Command=['C:\EnergyPlusV9-4-0\energyplus     -s C   -r   -c   -p EPNewFile        -w  EPWeather.epw      EPNewFile.idf'];
[status,cmdout]=system(Command); %Runs the EP simulation.

%Reads the raw simulator output from the output file (EPNewFile.csv).
ReadData=readmatrix(['EPNewFile.csv']);
SimTemperature=ReadData(:,[55 67]); %Raw simulator output (a matrix): Column i (where i = 1 or 2) gives the temperatures at the location of Sensor i at 24*Number_TimePoints_EveryHour equally-spaced time points, starting from 00:01 and ending at 24:00, over the time period 00:01 - 24:00 of a day.

%The following computes the average temperature for each of the two-hour intervals 00:01 - 02:00, 02:01 - 04:00, ..., 20:01 - 22:00 at each of the locations of Sensor 1 and Sensor 2.
Number_TimePoints_EveryTwoHours=Number_TimePoints_EveryHour*2;
OutputForCalibration=[];
for Sensor=1:2
    TemperatureOutput=SimTemperature(:,Sensor); 
    TemperatureOutput=reshape(TemperatureOutput,Number_TimePoints_EveryTwoHours,[]);
    AverageTemperature11Intervals=mean(TemperatureOutput); %Average temperature for each of 12 two-hour intervals, i.e., 00:01 - 02:00, 02:01 - 04:00, ..., 22:01 - 24:00.
    AverageTemperature11Intervals(end)=[]; %Removes the last value (average temperature for 22:01 to 24:00) to get 11 average temperatures.
    OutputForCalibration=[OutputForCalibration AverageTemperature11Intervals]; %Output vector, which has 2*11 elements.
end

end
%%
function Output=SimulatorPDE(x,Level)
MinInputs=[20  0.0143*0.8 0.2137*0.8];
MaxInputs=[100 0.0143*1.2 0.2137*1.2];
RangeInputs=MaxInputs-MinInputs;
Inputs=MinInputs+RangeInputs.*x;
ConvectionCoefficient=Inputs(1);
A=6.212;
B=Inputs(2);
D=363.56;
E=Inputs(3);
MassDensity=8170;

AmbientTemperature=1200; InitialTemperature=300;
ThermalConductivity=@(~,state) A+B*(state.u);              
SpecificHeat=@(~,state) D+E*(state.u);                                                   

W1=0.05*0.6;
H1=0.1*0.6;

thermalmodel=createpde('thermal','transient');
r1=[3 4 -W1 W1 W1 -W1 -H1 -H1 H1 H1];
gdm=r1';
gm=decsg(gdm,'R1',['R1']');
geometryFromEdges(thermalmodel,gm);

if Level==1
    mesh=generateMesh(thermalmodel,'Hmax',0.0125);
elseif Level==2
    mesh=generateMesh(thermalmodel,'Hmax',0.01);
end

for i=1:thermalmodel.Geometry.NumEdges
    nodes=findNodes(mesh,'region','Edge',i);
    if(all(abs(thermalmodel.Mesh.Nodes(1,nodes)+W1)<10^-6))
        LeftE=i;
    elseif(all(abs(thermalmodel.Mesh.Nodes(1,nodes)-W1)<10^-6))
        RightE=i;
    elseif(all(abs(thermalmodel.Mesh.Nodes(2,nodes)-H1)<10^-6))
        TopE=i;        
    elseif(all(abs(thermalmodel.Mesh.Nodes(2,nodes)+H1)<10^-6))
        BottomE=i;
    end
end
if(isempty(setdiff(1:thermalmodel.Geometry.NumEdges,[LeftE RightE TopE BottomE]))~=1)
    return
end
thermalProperties(thermalmodel,'ThermalConductivity',ThermalConductivity,'MassDensity',MassDensity,'SpecificHeat',SpecificHeat);
thermalBC(thermalmodel,'Edge',BottomE,'Temperature',InitialTemperature);
if Level==1
    thermalBC(thermalmodel,'Edge',[RightE LeftE TopE],'ConvectionCoefficient',ConvectionCoefficient,'AmbientTemperature',AmbientTemperature);
elseif Level==2
    thermalBC(thermalmodel,'Edge',[RightE LeftE TopE],'ConvectionCoefficient',ConvectionCoefficient,'AmbientTemperature',AmbientTemperature,'Emissivity',0.05);               
    thermalmodel.StefanBoltzmannConstant=5.670367e-8;                                                                                                                                                                                                 
end
thermalIC(thermalmodel,InitialTemperature);

tlist=0:90:900;
thermalresults=solve(thermalmodel,tlist);
nlevels=11;
standardizedlevels=linspace(-(nlevels-1)/nlevels,(nlevels-1)/nlevels,nlevels);
[xq,yq]=meshgrid(standardizedlevels*W1,standardizedlevels*H1);
Temperature=interpolateTemperature(thermalresults,xq,yq,2:numel(tlist));

Output=Temperature(:)';

end
%%
function Output=SimulatorEnv(x,Level)

MinInputs=[0.1 1];
MaxInputs=[2   34];
RangeInputs=MaxInputs-MinInputs;
Inputs=MinInputs+RangeInputs.*x; 
L=Inputs(1);
tau=Inputs(2);

D=0.04;
M=13;

s=1:0.4:3;
t=35:5:60;

ds=length(s);
dt=length(t);

ymatrix=zeros(dt,ds);

for ii=1:dt
    ti=t(ii);
    Ind=tau<ti;
    
    for jj=1:ds
        sj=s(jj);
        
        factor1a = M/sqrt(4*pi*D*ti);
        if Level==1
            factor1b = max(1+(-sj^2/(4*D*ti))/6,0)^6;
        elseif Level==2
            factor1b = exp(-sj^2/(4*D*ti));
        end
        term1 = factor1a*factor1b;
        
        term2 = 0;
        if Ind
            factor2a = M/sqrt(4*pi*D*(ti-tau));
            if Level==1
                factor2b = max(1+(-(sj-L)^2/(4*D*(ti-tau)))/6,0)^6;
            elseif Level==2
                factor2b = exp(-(sj-L)^2/(4*D*(ti-tau)));
            end
            term2 = factor2a*factor2b;
        end
        
        C = term1+term2;
        ymatrix(ii,jj) = sqrt(4*pi)*C;
    end
end

Output=ymatrix(:)';

end