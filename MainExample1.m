%Example 1
%Section 1: Sets input data for all model calibration methods.
%Dl is the initial LF design (an nl x d matrix) and Dh is the initial HF design (an nh x d matrix) for each bi-fidelity method in a trial;
%Yl is the initial LF output data (an nl x N matrix) and Yh is the initial HF output data (an nh x N matrix) for each bi-fidelity method in a trial; 
%Dh0 is the initial HF design (an nh0 x d matrix) for each single-fidelity method in a trial;
%Yh0 is the initial HF output data (an nh0 x N matrix) for each single-fidelity method in a trial;
%CostRatio is c_h/c_l; 
%Budget is the total budget for each bi-fidelity method in a trial; 
%Budget1 is the total budget for each single-fidelity method in a trial;
%w is the vector of field/physical data (a 1 x N vector).

%%Uncomment the commands below in this section if you want to generate a different set of 100 initial designs for each method, rather than use those 
%stored in the .mat file loaded by this script, to perform another set of 100 trials.
%{
clear all;clc,format compact 
SensorTemperature=[24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24.5;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24.5;24,24;24,24.5;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;23.5,24;23.5,24;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23;23.5,23;23,23;23,23;23,23;23,23;23,23;23,22.5;23,22.5;23,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22;22.5,22;22.5,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,21.5;22,21.5;22,21.5;22,21.5;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22.5;22.5,22;22.5,22;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22;22.5,22;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22,22.5;21,22.5;20,22.5;19,22.5;18,22.5;17.5,22;17,22;17,22;17.5,22;17.5,22;18,22;18,22;18.5,22;18.5,22;19,22;19,22;19.5,22;20,22;20,22;20.5,22;20.5,22;21,22;21,22;21.5,22;21.5,22;21.5,22;22,22;22,22;22,22;22,22;22,22.5;22,22.5;22,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;21.5,22.5;20.5,22.5;19,22.5;18.5,22.5;18,22.5;17.5,22;17,22;17.5,22;17.5,22;18,22;18,22;18.5,22;18.5,22;19,22;19,22;19.5,22;20,22;20,22;20.5,22;20.5,22;21,22;21,22;21.5,22;21.5,22;21.5,22;21.5,22;22,22;22,22;22,22;22,22;22,22.5;22,22.5;22,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22,22.5;21,22.5;20,22.5;19,22.5;18,22.5;17.5,22;17,22;17,22;17.5,22;17.5,22;18,22;18,22;18.5,22;18.5,22;19,22;19,22;19.5,22;20,22;20,22;20.5,22;20.5,22;21,22;21,22;21,22;21.5,22;21.5,22;21.5,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;21.5,22.5;20.5,22.5;19.5,22.5;18.5,22.5;18,22.5;17.5,22;17,22;17,22;17.5,22;17.5,22;18,22;18,22;18.5,22;18.5,22;19,22;19,22;19.5,22;20,22;20,22;20.5,22;20.5,22;21,22;21,22;21,22;21.5,22;21.5,22;21.5,22;21.5,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;23,23;23,23;23,23;23,23;23,23;23,23;23,23;23,23;23.5,23;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5];
SensorTemperatureEveryTwoHours=reshape(SensorTemperature,120,[]);
w=mean(SensorTemperatureEveryTwoHours);

d=3;
Case=1;
xstar='unknown';
nl=16;
nh=8;
nh0=12;
CostRatio=4;
InitialBudget=nl*1+CostRatio*nh;
InitialBudget0=nh0*CostRatio;
if(InitialBudget~=InitialBudget0)
    return
end
Budget=InitialBudget0+15;
Budget1=Budget+1;

parfor id=1:100
    id
    [Dl,Dh]=GenerateNestedLHD(nl,nh,d,1e5);     
    [Dh0]=GenerateNestedLHD(nh0,nh0,d,1e5);     
    
    Dls(:,:,id)=Dl;
    Dhs(:,:,id)=Dh;
    Dh0s(:,:,id)=Dh0;    
end

for id=1:100
    id
    Dl=Dls(:,:,id);
    Dh=Dhs(:,:,id);
    Dh0=Dh0s(:,:,id);
    clear Yl Yh
    for jd=1:nl
        Yl(jd,:)=Simulator(Dl(jd,:),1,Case);
    end
    for jd=1:nh
        Yh(jd,:)=Simulator(Dh(jd,:),2,Case);
    end
    clear Yh0
    for jd=1:nh0
        Yh0(jd,:)=Simulator(Dh0(jd,:),2,Case);
    end
    
    MultiDataInput(id).Dl=Dl;               MultiDataInput(id).Yl=Yl;
    MultiDataInput(id).Dh=Dh;               MultiDataInput(id).Yh=Yh;
    MultiDataInput(id).xstar=xstar;
    MultiDataInput(id).w=w;                 MultiDataInput(id).CostRatio=CostRatio;
    MultiDataInput(id).Budget=Budget;       MultiDataInput(id).Case=Case;
    
    SingleDataInput(id).Dl=[];              SingleDataInput(id).Yl=[];
    SingleDataInput(id).Dh=Dh0;             SingleDataInput(id).Yh=Yh0;
    SingleDataInput(id).xstar=xstar;
    SingleDataInput(id).w=w;                SingleDataInput(id).CostRatio=CostRatio;
    SingleDataInput(id).Budget=Budget1;     SingleDataInput(id).Case=Case;
    
end
%}
%save Example1InputData.mat
%clear all
%load Example1InputData.mat MultiDataInput SingleDataInput
%save Example1InputData.mat
clear all;clc,format compact 
%Column i (i=1,2) of SensorTemperature gives the temperature measurements made by Sensor i every minute from 12:01 a.m. to 10:00 p.m. of a day.
SensorTemperature=[24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24.5;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24.5;24,24;24,24.5;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;23.5,24;23.5,24;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23;23.5,23;23,23;23,23;23,23;23,23;23,23;23,22.5;23,22.5;23,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22;22.5,22;22.5,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,21.5;22,21.5;22,21.5;22,21.5;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22.5;22.5,22;22.5,22;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22;22.5,22;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22,22.5;21,22.5;20,22.5;19,22.5;18,22.5;17.5,22;17,22;17,22;17.5,22;17.5,22;18,22;18,22;18.5,22;18.5,22;19,22;19,22;19.5,22;20,22;20,22;20.5,22;20.5,22;21,22;21,22;21.5,22;21.5,22;21.5,22;22,22;22,22;22,22;22,22;22,22.5;22,22.5;22,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;21.5,22.5;20.5,22.5;19,22.5;18.5,22.5;18,22.5;17.5,22;17,22;17.5,22;17.5,22;18,22;18,22;18.5,22;18.5,22;19,22;19,22;19.5,22;20,22;20,22;20.5,22;20.5,22;21,22;21,22;21.5,22;21.5,22;21.5,22;21.5,22;22,22;22,22;22,22;22,22;22,22.5;22,22.5;22,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22,22.5;21,22.5;20,22.5;19,22.5;18,22.5;17.5,22;17,22;17,22;17.5,22;17.5,22;18,22;18,22;18.5,22;18.5,22;19,22;19,22;19.5,22;20,22;20,22;20.5,22;20.5,22;21,22;21,22;21,22;21.5,22;21.5,22;21.5,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;21.5,22.5;20.5,22.5;19.5,22.5;18.5,22.5;18,22.5;17.5,22;17,22;17,22;17.5,22;17.5,22;18,22;18,22;18.5,22;18.5,22;19,22;19,22;19.5,22;20,22;20,22;20.5,22;20.5,22;21,22;21,22;21,22;21.5,22;21.5,22;21.5,22;21.5,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;23,23;23,23;23,23;23,23;23,23;23,23;23,23;23,23;23.5,23;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5];
SensorTemperatureEveryTwoHours=reshape(SensorTemperature,120,[]);
w=mean(SensorTemperatureEveryTwoHours);

d=3;
Case=1;
xstar='unknown';
nl=16;
nh=8;
nh0=12;
CostRatio=4;
InitialBudget=nl*1+CostRatio*nh;
InitialBudget0=nh0*CostRatio;
if(InitialBudget~=InitialBudget0)
    return
end
Budget=InitialBudget0+15;
Budget1=Budget+1;

load Example1InputData.mat MultiDataInput SingleDataInput
%%
%Section 2: Runs all model calibration methods.
Z_BC=1;    Z_ID=0;   Z_SR=2;
ZMLFSSE=1; ZLFSSE=0; AccuracyLevel=1; t=0.99; 
for id=1:100
    id
    [T_MBC_AGP{id,1},Data_MBC_AGP{id,1},RunTime_MBC_AGP(id)]=CalibrationAGP(MultiDataInput(id),Z_BC,ZMLFSSE,AccuracyLevel); 'MBC-AGP'
    [T_BC_AGP{id,1},Data_BC_AGP{id,1},RunTime_BC_AGP(id)]=CalibrationAGP(MultiDataInput(id),Z_BC,ZLFSSE,AccuracyLevel); 'BC-AGP'
    [T_MID_AGP{id,1},Data_MID_AGP{id,1},RunTime_MID_AGP(id)]=CalibrationAGP(MultiDataInput(id),Z_ID,ZMLFSSE,AccuracyLevel); 'MID-AGP'
    [T_SR_AGP{id,1},Data_SR_AGP{id,1},RunTime_SR_AGP(id)]=CalibrationAGP(MultiDataInput(id),Z_SR,ZLFSSE,AccuracyLevel); 'SR-AGP'
    [T_Nested{id,1},Data_Nested{id,1},RunTime_Nested(id)]=CalibrationNested(MultiDataInput(id),AccuracyLevel); 'Nested'
    [T_SVD_AGP{id,1},Data_SVD_AGP{id,1},RunTime_SVD_AGP(id)]=CalibrationSVDAGP(MultiDataInput(id),AccuracyLevel,t); 'SVD-AGP'
    [T_BC_GP{id,1},Data_BC_GP{id,1},RunTime_BC_GP(id)]=CalibrationBCGP(SingleDataInput(id),AccuracyLevel); 'BC-GP'
    [T_SR_GP{id,1},Data_SR_GP{id,1},RunTime_SR_GP(id)]=CalibrationSRGP(SingleDataInput(id),AccuracyLevel); 'SR-GP'
    [T_SVD{id,1},Data_SVD{id,1},RunTime_SVD(id)]=CalibrationSVD(SingleDataInput(id),AccuracyLevel,t); 'SVD'    
    save Example1.mat
end
%%
%Section 3: Constructs figures and a table to illustrate results.
clear all;clc,format compact 
load Example1.mat
idx=(1:100);
Labels={'MBC-AGP','BC-AGP','MID-AGP','SR-AGP','Nested','SVD-AGP','BC-GP','SR-GP','SVD'}';
RecordTable=[T_MBC_AGP(idx) T_BC_AGP(idx) T_MID_AGP(idx) T_SR_AGP(idx) T_Nested(idx) T_SVD_AGP(idx) T_BC_GP(idx) T_SR_GP(idx) T_SVD(idx)];
RecordData=[Data_MBC_AGP(idx) Data_BC_AGP(idx) Data_MID_AGP(idx) Data_SR_AGP(idx) Data_Nested(idx) Data_SVD_AGP(idx) Data_BC_GP(idx) Data_SR_GP(idx) Data_SVD(idx)];

xstarML=[0.458055953979492 0.971427917480469 0.999999923706055];
ShxstarML=3.52649258492504;

MinInputs=[0.12 0.0004 0];
MaxInputs=[0.3 0.001 0.975];
RangeInputs=MaxInputs-MinInputs;
UnnormalizedxstarML=MinInputs+RangeInputs.*xstarML;

for idxMethod=9:-1:1

    for idxTrial=1:numel(idx)
        Table=RecordTable{idxTrial,idxMethod};
        
        if idxMethod<=2 || idxMethod==7
            phiEnd(idxTrial,idxMethod)=Table.phis(end,:);
        end
        
        costs=[1 CostRatio]';
        ShxhatstarMLs=Table.ShxhatstarMLs;
        xhatstarMLs=Table.xhatstarMLs;

        L2xhatstarMLs=sum((xhatstarMLs-xstarML).^2,2).^0.5;
        Levels=Table.Level;
        Costs=cumsum(costs(Levels));

        ShxhatstarMLsEnd(idxTrial,idxMethod)=ShxhatstarMLs(end,:);
        L2xhatstarMLsEnd(idxTrial,idxMethod)=L2xhatstarMLs(end);
        
        if idxMethod==7 || idxMethod==8 || idxMethod==9
            InterpolatedShxhatstarMLs(1:Budget1,idxMethod,idxTrial)=interp1(Costs,ShxhatstarMLs,1:Budget1);
            InterpolatedL2xhatstarMLs(1:Budget1,idxMethod,idxTrial)=interp1(Costs,L2xhatstarMLs,1:Budget1);
        elseif idxMethod<=4
            InterpolatedShxhatstarMLs(1:Budget,idxMethod,idxTrial)=interp1(Costs,ShxhatstarMLs,1:Budget);
            InterpolatedL2xhatstarMLs(1:Budget,idxMethod,idxTrial)=interp1(Costs,L2xhatstarMLs,1:Budget);
        elseif idxMethod==5 || idxMethod==6
            deleteidx=(nl+nh+1):2:size(Table,1);
            Costs(deleteidx,:)=[];
            ShxhatstarMLs(deleteidx,:)=[];
            L2xhatstarMLs(deleteidx,:)=[];
            
            InterpolatedShxhatstarMLs(1:Budget,idxMethod,idxTrial)=interp1(Costs,ShxhatstarMLs,1:Budget);
            InterpolatedL2xhatstarMLs(1:Budget,idxMethod,idxTrial)=interp1(Costs,L2xhatstarMLs,1:Budget);
            
        end
    end
end
meanInterpolatedShxhatstarMLs=mean(InterpolatedShxhatstarMLs,3);
meanInterpolatedShxhatstarMLsminusShxstarML=meanInterpolatedShxhatstarMLs-ShxstarML;
meanInterpolatedL2xhatstarMLs=mean(InterpolatedL2xhatstarMLs,3);

idx1=1;
for idx2=1:9
    [~,ttest_pval_Sh(idx2,1)]=ttest(ShxhatstarMLsEnd(:,idx1),ShxhatstarMLsEnd(:,idx2));
    [~,ttest_pval_L2(idx2,1)]=ttest(L2xhatstarMLsEnd(:,idx1),L2xhatstarMLsEnd(:,idx2));
end
AverageSh=mean(ShxhatstarMLsEnd)'; AverageL2=mean(L2xhatstarMLsEnd)';
disp('Table H.1')
Table1=table(Labels,AverageSh,ttest_pval_Sh,AverageL2,ttest_pval_L2)

htmlGray=[128 128 128]/255;
htmlGreen=[0.4660 0.6740 0.1880];

Jump=5;
for i=1:9
    if(i<=6)
        JJ(i)= mean(log(meanInterpolatedShxhatstarMLsminusShxstarML(InitialBudget:Budget,i)));
    else
        JJ(i)= mean(log(meanInterpolatedShxhatstarMLsminusShxstarML(InitialBudget:Budget1,i)));        
    end
end
[~,indicesJJ]=sort(JJ);
for i=1:9
    Shift=mod(find(indicesJJ==i),Jump);
    if(i<=6)
        II{i}=unique([InitialBudget (InitialBudget+Shift):Jump:Budget Budget]);
    else
        II{i}=unique([InitialBudget (InitialBudget+Shift):Jump:Budget1 Budget1]);        
    end
end

figure,clf
tiledlayout(1,2,'Padding','none','TileSpacing','none');
nexttile
FontSize1=24;
linewidth=4;
MarkerSize1=15;
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,1),'ko-','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{1}),hold on
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,2),'b:o','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerFaceColor','b','MarkerIndices',II{2})
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,3),'k^-','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{3})
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,4),'--v','linewidth',linewidth,'color',htmlGray,'MarkerSize',MarkerSize1,'MarkerIndices',II{4})
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,5),':s','linewidth',linewidth,'color',htmlGreen,'MarkerFaceColor',htmlGreen,'MarkerSize',MarkerSize1,'MarkerIndices',II{5})
plot(1:Budget,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget,6),'b-x','linewidth',linewidth,'MarkerSize',MarkerSize1+10,'MarkerIndices',II{6})
plot(1:Budget1,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget1,7),':rs','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{7})
plot(1:Budget1,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget1,8),'--h','linewidth',linewidth,'color',[0.00,0.45,0.74],'MarkerFaceColor',[0.00,0.45,0.74],'MarkerSize',MarkerSize1,'MarkerIndices',II{8})
plot(1:Budget1,meanInterpolatedShxhatstarMLsminusShxstarML(1:Budget1,9),':d','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{9})
xlabel('Computational cost');
set(gca,'YScale','log','FontSize',FontSize1,'FontWeight','Bold','LineWidth',3);
ylabel('Average $S_h(\hat{\textbf{x}}^*_{\mathbf{ML}})-$3.526493','Interpreter','latex','FontSize',32);
leg=legend(Labels,'NumColumns',3,'Location','northeast');
leg.ItemTokenSize=[74,50];
yticks([0.0625 0.125 0.25 0.5 1 2 4 8])
ylim([0.05 12.2])
xticks([InitialBudget:2:Budget1])
xlim([InitialBudget-0.1 Budget1+0.1])
title('(a)','FontWeight','Bold')

for i=1:9
    if(i<=6)
        JJ(i)= mean(meanInterpolatedL2xhatstarMLs(InitialBudget:Budget,i));
    else
        JJ(i)= mean(meanInterpolatedL2xhatstarMLs(InitialBudget:Budget1,i));        
    end
end
[~,indicesJJ]=sort(JJ);
for i=1:9
    Shift=mod(find(indicesJJ==i),Jump);
    if(i<=6)
        II{i}=unique([InitialBudget (InitialBudget+Shift):Jump:Budget Budget]);
    else
        II{i}=unique([InitialBudget (InitialBudget+Shift):Jump:Budget1 Budget1]);        
    end
end

nexttile
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,1),'ko-','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{1}),hold on
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,2),'b:o','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerFaceColor','b','MarkerIndices',II{2})
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,3),'k^-','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{3})
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,4),'--v','linewidth',linewidth,'color',htmlGray,'MarkerSize',MarkerSize1,'MarkerIndices',II{4})
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,5),':s','linewidth',linewidth,'color',htmlGreen,'MarkerFaceColor',htmlGreen,'MarkerSize',MarkerSize1,'MarkerIndices',II{5})
plot(1:Budget,meanInterpolatedL2xhatstarMLs(1:Budget,6),'b-x','linewidth',linewidth,'MarkerSize',MarkerSize1+10,'MarkerIndices',II{6})
plot(1:Budget1,meanInterpolatedL2xhatstarMLs(1:Budget1,7),':rs','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{7})
plot(1:Budget1,meanInterpolatedL2xhatstarMLs(1:Budget1,8),'--h','linewidth',linewidth,'color',[0.00,0.45,0.74],'MarkerFaceColor',[0.00,0.45,0.74],'MarkerSize',MarkerSize1,'MarkerIndices',II{8})
plot(1:Budget1,meanInterpolatedL2xhatstarMLs(1:Budget1,9),':d','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{9})
xlabel('Computational cost');
set(gca,'FontWeight','bold','FontSize',FontSize1);
ylabel('Average $L_2(\hat{\textbf{x}}^*_{\mathbf{ML}})$','Interpreter','latex','FontSize',32);
leg=legend(Labels,'NumColumns',3,'Location','northeast');
leg.ItemTokenSize=[74,50];
set(findobj(gcf,'type','axes'),'FontWeight','Bold','LineWidth',3);
yticks([0.2:0.1:0.8])
ylim([0.19 0.825])
xticks([InitialBudget:2:Budget1])
xlim([InitialBudget-0.1 Budget1+0.1])
title('(b)','FontWeight','Bold')
set(gcf,'Position',[0 100 1920 615])

figure,clf
Labels2Method={'MBC-AGP','BC-AGP','BC-GP'};
boxplot(phiEnd(:,[1 2 7]),'Labels',Labels2Method,'OutlierSize',10,'Widths',0.8*[1 1 1])
set(findobj(gca,'type','line'),'LineWidth',2)
set(findobj(gcf,'type','axes'),'FontSize',27,'FontWeight','Bold','LineWidth',3);
ylabel('$ \hat \varphi$','Interpreter','latex','FontSize',50,'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','baseline')
set(gca,'Position',[0.15 0.15 0.83 0.83])
yticks([-0.1:0.1:0.6])
set(gca,'yGrid','on','GridLineStyle','--')
ylim([-0.1018 0.6365])
set(gcf,'Position',[109 159 900 372])
medians=median(phiEnd(:,[1 2 7]));
FontSize2=21;
text(1,1.085*medians(1),['Median=' num2str(medians(1),2)],'HorizontalAlignment','center','FontSize',FontSize2,'FontWeight','Bold')
text(2,1.08*medians(2),['Median=' num2str(medians(2),2)],'HorizontalAlignment','center','FontSize',FontSize2,'FontWeight','Bold')
text(3,1.15*medians(3),['Median=' num2str(medians(3),2)],'HorizontalAlignment','center','FontSize',FontSize2,'FontWeight','Bold')
xlim([0.45 3.55])

idxTrial=54;
pairindices=[1 2; 1 3; 2 3];
for idxMethod=1:2
    Data=RecordData{idxTrial,idxMethod};
    Table=RecordTable{idxTrial,idxMethod};
    xhatstarMLsEnd=Table.xhatstarMLs(end,:);

    InitialDh=Data.Dh(1:nh,:);
    InitialDl=Data.Dl(1:nl,:);
    FollowDh=Data.Dh(nh+1:end,:);
    FollowDl=Data.Dl(nl+1:end,:);
    
    figure,clf
    tiledlayout(1,3,'Padding','none','TileSpacing','none');
    
    for kd=1:3
        pd1=pairindices(kd,1);
        pd2=pairindices(kd,2);
        nexttile
        Markersize1=15;
        linewidth=1.5;
        plot(InitialDh(:,pd1),InitialDh(:,pd2),'bs','linewidth',linewidth,'markersize',Markersize1)
        hold on
        plot(InitialDl(:,pd1),InitialDl(:,pd2),'bx','linewidth',linewidth,'markersize',Markersize1)
        
        plot(FollowDh(:,pd1),FollowDh(:,pd2),'ko','linewidth',linewidth,'markersize',Markersize1)
        
        plot(FollowDl(:,pd1),FollowDl(:,pd2),'k+','linewidth',linewidth,'markersize',Markersize1)
        
        plot(xstarML(:,pd1),xstarML(:,pd2),'kp','MarkerSize',25,'linewidth',linewidth)
        plot(xhatstarMLsEnd(:,pd1),xhatstarMLsEnd(:,pd2),'k^','MarkerSize',25,'linewidth',linewidth)
        
        xlabel(['x_' num2str(pd1)])
        ylabel(['x_' num2str(pd2)],'Rotation',0,'HorizontalAlignment','right')

        xticks([0:0.2:1])
        yticks([0:0.2:1])

        lim0=0.02;
        xlim([-lim0 1+lim0])
        ylim([-lim0 1+lim0])
                         
    end
    
    set(findobj(gcf,'type','axes'),'FontSize',17,'FontWeight','Bold','LineWidth',1);
    if idxMethod==1
        sgtitle('(a) MBC-AGP','fontsize',25,'FontWeight','Bold')
        set(gcf,'Position',[0 450 1600 350])        
    else
        sgtitle('(b) BC-AGP','fontsize',25,'FontWeight','Bold')
        set(gcf,'Position',[0 0 1600 350])        
    end
    
end

idxTrial=47;
for idxMethod=1:9
    Table=RecordTable{idxTrial,idxMethod};
    yhxhatstarMLs(idxMethod,:)=RecordData{idxTrial,idxMethod}.yhxhatstarMLs(end,:);
    ShxhatstarMLsEndH2H3(idxMethod,1)=Table.ShxhatstarMLs(end,:);
end

linewidth=3;
MarkerSize1=15;
Labels1={'MBC-AGP','BC-AGP','MID-AGP','SR-AGP','Nested','SVD-AGP','Field data'};
html1=[0.3010 0.7450 0.9330];

II=[]; 
for kd=1:2    
    kdidx=(kd-1)*11+[1:11];    
    JJ=[];
for i=2:6
    JJ(i-1)= mean(yhxhatstarMLs(i,kdidx));
end
[~,indicesJJ]=sort(JJ);
for i=2:6
    Jump=2;
    Shift=mod(find(indicesJJ==(i-1)),Jump);
    II{i,kd}=[(1+Shift):Jump:11];
end
    JJ2=[];
for i=7:9
    JJ2(i-6)=mean(yhxhatstarMLs(i,kdidx));
end
[~,indicesJJ2]=sort(JJ2);
for i=7:9
    Jump=2;
    Shift=mod(find(indicesJJ2==(i-6)),Jump);
    II{i,kd}=[(1+Shift):Jump:11];
end
    II{10,kd}=1:11;
end

figure,clf
for kd=1:2
    subplot(1,2,kd)
    
    kdidx=(kd-1)*11+[1:11];
    plot(yhxhatstarMLs(1,kdidx),'k:','linewidth',linewidth+3,'MarkerSize',MarkerSize1,'MarkerIndices',II{1,kd})
    hold on
    plot(yhxhatstarMLs(2,kdidx),':v','color','m','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{2,kd})
    plot(yhxhatstarMLs(3,kdidx),':ro','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{3,kd})
    plot(yhxhatstarMLs(4,kdidx),'--bs','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{4,kd})
    plot(yhxhatstarMLs(5,kdidx),'--d','color',html1,'linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{5,kd})
    plot(yhxhatstarMLs(6,kdidx),'--r^','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{6,kd})
    plot(w(kdidx),'-p','color',htmlGreen,'linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{10,kd})
    xticks(1:11)
    
    if kd==1
        ylim([20.8 24.6])
        title('(a)','FontWeight','bold')
        set(gca,'Position',[0.065 0.215 0.425 0.745])
    else
        ylim([21.4 24.6])
        title('(b)','FontWeight','bold')
        set(gca,'Position',[0.56 0.215 0.425 0.745])
    end
    
    set(findobj(gca,'type','axes'),'FontWeight','Bold','LineWidth',2);
    bp=gca; bp.FontSize=20;
    bp.XAxis.FontSize=20;
    bp.YAxis.FontSize=23;
    
    row1={'12:01 a.m.','2:01 a.m.','4:01 a.m.','6:01 a.m.','  8:01 a.m.','  10:01 a.m.','12:01 p.m.','2:01 p.m.','4:01 p.m.','6:01 p.m.','  8:01 p.m.'};
    row2={'to 2:00 a.m.','to 4:00 a.m.','to 6:00 a.m.','to 8:00 a.m.','to 10:00 a.m.','to 12:00 p.m.','to 2:00 p.m.','to 4:00 p.m.','to 6:00 p.m.','to 8:00 p.m.','to 10:00 p.m.'};
    labelArray=[row1;row2];
    tickLabels=strtrim(sprintf('%s\\newline%s\n',labelArray{:}));
    bp.XTickLabel=tickLabels;
    bp.XTickLabelRotation=90;
    
    xlabel('Time interval')
    
    xlim([0.8 11.2])
    ylabel('Temperature (Celsius)')
    leg=legend(Labels1,'NumColumns',2);
    leg.ItemTokenSize=[60,50];
    set(gcf,'Position',[0 59 1800 963])
end

Labels2={'MBC-AGP','BC-GP','SR-GP','SVD','Field data'};

figure,clf
for kd=1:2
    subplot(1,2,kd)
    
    kdidx=(kd-1)*11+[1:11];
    plot(yhxhatstarMLs(1,kdidx),'k:','linewidth',linewidth+3,'MarkerSize',MarkerSize1,'MarkerIndices',II{1,kd})
    hold on
    plot(yhxhatstarMLs(7,kdidx),'<--b','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{7,kd})
    plot(yhxhatstarMLs(8,kdidx),':r>','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{8,kd})
    plot(yhxhatstarMLs(9,kdidx),'-.ks','linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{9,kd})
    plot(w(kdidx),'-p','color',htmlGreen,'linewidth',linewidth,'MarkerSize',MarkerSize1,'MarkerIndices',II{10,kd})
    xticks(1:11)
    
    if kd==1
        ylim([20.8 24.6])
        title('(a)','FontWeight','bold')
        set(gca,'Position',[0.065 0.255 0.425 0.695])
    else
        ylim([21.4 24.6])
        title('(b)','FontWeight','bold')
        set(gca,'Position',[0.56 0.255 0.425 0.695])       
    end
    
    set(findobj(gca,'type','axes'),'FontWeight','Bold','LineWidth',2);
    bp=gca; bp.FontSize=20;
    bp.XAxis.FontSize=20;
    bp.YAxis.FontSize=23;
    
    row1={'12:01 a.m.','2:01 a.m.','4:01 a.m.','6:01 a.m.','  8:01 a.m.','  10:01 a.m.','12:01 p.m.','2:01 p.m.','4:01 p.m.','6:01 p.m.','  8:01 p.m.'};
    row2={'to 2:00 a.m.','to 4:00 a.m.','to 6:00 a.m.','to 8:00 a.m.','to 10:00 a.m.','to 12:00 p.m.','to 2:00 p.m.','to 4:00 p.m.','to 6:00 p.m.','to 8:00 p.m.','to 10:00 p.m.'};
    labelArray=[row1;row2];
    tickLabels=strtrim(sprintf('%s\\newline%s\n',labelArray{:}));
    bp.XTickLabel=tickLabels;
    bp.XTickLabelRotation=90;
    
    xlabel('Time interval')
    
    xlim([0.8 11.2])
    ylabel('Temperature (Celsius)')
    leg=legend(Labels2,'NumColumns',2);
    leg.ItemTokenSize=[60,50];
    set(gcf,'Position',[0 59 1800 800])
end
%%
%This section gives results for Example 1 reported in the paragraph containing equation (H.1) in Appendix H, which are found by computing various integrals based on a quadrature rule.
clear all;clc,format compact;format long g 
load Example1GridData26.mat  AllYh AllYl w GridPoints
load Example1.mat MultiDataInput
AllSh=sum((AllYh-w).^2,2);
AllSl=sum((AllYl-w).^2,2);

phi=0.38;
% GridPoints: 26^3 by 3 matrix.
% AllYh: 26^3 by 22 matrix.
% AllYl: 26^3 by 22 matrix.
nlevels=26;
nlevelsm1=nlevels-1;
nlevelsm2=nlevels-2;
cellno=1;
for a=0:(1/nlevelsm1):nlevelsm2*(1/nlevelsm1)
    for b=0:(1/nlevelsm1):nlevelsm2*(1/nlevelsm1)
        for c=0:(1/nlevelsm1):nlevelsm2*(1/nlevelsm1)
            Lowervertex=[a,b,c];
            V=(fullfact([2 2 2])-1)*(1/nlevelsm1);
            Allvertices{cellno}=repmat(Lowervertex,8,1)+V;
            cellno=cellno+1;
        end
    end
end

for k=1:nlevelsm1^3
    [mindists,minidx]=min(pdist2(Allvertices{k},GridPoints),[],2);
    Cellverticesindices(k,:)=minidx;
    CellPoints=GridPoints(Cellverticesindices(k,:),:);
    vol=prod(range(CellPoints,1));
    if((size(CellPoints,1)~=8) || (abs(vol-(1/nlevelsm1)^3)>10^-12) || (min(pdist(CellPoints))<10^-12) || (max(mindists)>10^-12))
        input('error')
        return
    end
end

AllSh_Vertices=AllSh(Cellverticesindices);
AllSl_Vertices=AllSl(Cellverticesindices);

Sh2_Vertices=AllSh_Vertices.^2;
Mean_Sh2=mean(Sh2_Vertices,'all');

L2_ShSl_Vertices=(AllSh_Vertices-AllSl_Vertices).^2;
L2_ShSl=mean(L2_ShSl_Vertices,'all');

ID_BC_or_SR=1;
Zh_Vertices=TransformData(AllSh_Vertices,phi,ID_BC_or_SR);
Zl_Vertices=TransformData(AllSl_Vertices,phi,ID_BC_or_SR);

ZhZl_Vertices=Zh_Vertices.*Zl_Vertices;
Mean_ZhZl=mean(ZhZl_Vertices,'all');

Mean_Zh=mean(Zh_Vertices,'all');
Mean_Zl=mean(Zl_Vertices,'all');

Zh2_Vertices=Zh_Vertices.^2;
Mean_Zh2=mean(Zh2_Vertices,'all');
Var_Zh=Mean_Zh2-(Mean_Zh)^2;

Zl2_Vertices=Zl_Vertices.^2;
Mean_Zl2=mean(Zl2_Vertices,'all');
Var_Zl=Mean_Zl2-(Mean_Zl)^2;

NoTrials=100;
L2_ShSlPlus=zeros(NoTrials,1); 
corr_ZhZlPlus=zeros(NoTrials,1); 
NormalizedL2_ShSlPlus=zeros(NoTrials,1);
for Trial=1:NoTrials
    Dl=MultiDataInput(Trial).Dl; Dh=MultiDataInput(Trial).Dh;
    Yl=MultiDataInput(Trial).Yl; Yh=MultiDataInput(Trial).Yh;
    [AllYlModified,ahati_bhati]=regress_aibi(Dl,Dh,Yl,Yh,AllYl);
    AllSlPlus=sum((AllYlModified-w).^2,2);
    
    AllSlPlus_Vertices=AllSlPlus(Cellverticesindices);
    L2_ShSlPlus_Vertices=(AllSh_Vertices-AllSlPlus_Vertices).^2;
    L2_ShSlPlus(Trial,1)=mean(L2_ShSlPlus_Vertices,'all');
    
    ZlPlus_Vertices=TransformData(AllSlPlus_Vertices,phi,ID_BC_or_SR);
    
    ZhZlPlus_Vertices=Zh_Vertices.*ZlPlus_Vertices;
    Mean_ZhZlPlus=mean(ZhZlPlus_Vertices,'all');
    
    Mean_ZlPlus=mean(ZlPlus_Vertices,'all');
    
    ZlPlus2_Vertices=ZlPlus_Vertices.^2;
    Mean_ZlPlus2=mean(ZlPlus2_Vertices,'all');
    Var_ZlPlus=Mean_ZlPlus2-(Mean_ZlPlus)^2;
    
    corr_ZhZlPlus(Trial,1)=(Mean_ZhZlPlus-Mean_Zh*Mean_ZlPlus)/(Var_Zh^0.5*Var_ZlPlus^0.5);
    
    NormalizedL2_ShSlPlus(Trial,1)=(L2_ShSlPlus(Trial,1)/Mean_Sh2)^0.5;
end

%%%%%%Results in Appendix H

%‖S_h (∙)-S_l (∙)‖_2/‖S_h (∙)‖_2.
disp('‖S_h (∙)-S_l (∙)‖_2/‖S_h (∙)‖_2 :')
NormalizedL2_ShSl=(L2_ShSl/Mean_Sh2)^0.5

%The 0.5 and 0.95 quantiles of ‖S_h (∙)-S_l^+ (∙)‖_2/‖S_h (∙)‖_2.
disp('The 0.5 and 0.95 quantiles of ‖S_h (∙)-S_l^+ (∙)‖_2/‖S_h (∙)‖_2 :')
prctile(NormalizedL2_ShSlPlus,[50 95])

%0.05 quantile and median of the correlation between g_φ (S_l^+ (X)) and g_φ (S_h (X)) over the 100 trials.
disp('0.05 quantile and median of correlation between g_φ (S_l^+ (X)) and g_φ (S_h (X)) :')
TwoQuantilesofcorr_ZhZlPlus=prctile(corr_ZhZlPlus,[5 50]) 

%The correlation between g_φ (S_h (X)) and g_φ (S_l (X)).
disp('Correlation between g_φ (S_h (X)) and g_φ (S_l (X)) :')
corr_ZhZl=(Mean_ZhZl-Mean_Zh*Mean_Zl)/(Var_Zh^0.5*Var_Zl^0.5)

disp('Result of Equation (H.1) :')
H1=(1-TwoQuantilesofcorr_ZhZlPlus(2)^2)/(1-corr_ZhZl^2)

function [YlTestModified,ahati_bhati]=regress_aibi(Dl,Dh,Yl,Yh,YlTest)
N=size(Yl,2);
[~,idxinDl,idxinDh]=intersect(Dl,Dh,'rows','stable');
SameYl=Yl(idxinDl,:);
SameYh=Yh(idxinDh,:);
OnesVec=ones(numel(idxinDh),1);
ahati_bhati=zeros(2,N);
for kd=1:N
    ModelMatrix=[OnesVec,SameYl(:,kd)];
    lastwarn('');
    ahati_bhati(:,kd)=regress(SameYh(:,kd),ModelMatrix);
    [warnMsg,~]=lastwarn;
    if contains(warnMsg,'X is rank deficient to within machine precision.')
        return
    end
end
YlTestModified=ahati_bhati(1,:)+YlTest.*ahati_bhati(2,:);
end