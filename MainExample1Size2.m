%Section 1: Sets parameters for all calibration methods 
clear,clc,format compact 
%Temperature measurement over 12:01 a.m. to 10:00 p.m.
SensorTemperature=[24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24.5,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24.5;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24.5;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24.5;24,24;24,24.5;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;24,24;23.5,24;23.5,24;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23;23.5,23;23,23;23,23;23,23;23,23;23,23;23,22.5;23,22.5;23,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22;22.5,22;22.5,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;21.5,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,21.5;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,21.5;22,21.5;22,21.5;22,21.5;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22.5;22.5,22;22.5,22;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22;22.5,22;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22,22.5;21,22.5;20,22.5;19,22.5;18,22.5;17.5,22;17,22;17,22;17.5,22;17.5,22;18,22;18,22;18.5,22;18.5,22;19,22;19,22;19.5,22;20,22;20,22;20.5,22;20.5,22;21,22;21,22;21.5,22;21.5,22;21.5,22;22,22;22,22;22,22;22,22;22,22.5;22,22.5;22,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;21.5,22.5;20.5,22.5;19,22.5;18.5,22.5;18,22.5;17.5,22;17,22;17.5,22;17.5,22;18,22;18,22;18.5,22;18.5,22;19,22;19,22;19.5,22;20,22;20,22;20.5,22;20.5,22;21,22;21,22;21.5,22;21.5,22;21.5,22;21.5,22;22,22;22,22;22,22;22,22;22,22.5;22,22.5;22,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22,22.5;21,22.5;20,22.5;19,22.5;18,22.5;17.5,22;17,22;17,22;17.5,22;17.5,22;18,22;18,22;18.5,22;18.5,22;19,22;19,22;19.5,22;20,22;20,22;20.5,22;20.5,22;21,22;21,22;21,22;21.5,22;21.5,22;21.5,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;21.5,22.5;20.5,22.5;19.5,22.5;18.5,22.5;18,22.5;17.5,22;17,22;17,22;17.5,22;17.5,22;18,22;18,22;18.5,22;18.5,22;19,22;19,22;19.5,22;20,22;20,22;20.5,22;20.5,22;21,22;21,22;21,22;21.5,22;21.5,22;21.5,22;21.5,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22.5,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22;22,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;22.5,22.5;23,23;23,23;23,23;23,23;23,23;23,23;23,23;23,23;23.5,23;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5;23.5,23.5];

Dim=3;
XTrue='unkown';
nh=12;
nl=24;
nh0=18;
RatioCost=4;
Case=1;
InitialBudget = nl*1+RatioCost*nh
InitialBudget0= nh0*RatioCost
Budget=InitialBudget0+15
Budget1=Budget+1;

SensorTemperatureEveryTwoHours=reshape(SensorTemperature,120,[]);
PhysData=mean(SensorTemperatureEveryTwoHours);
load Example1Size2.mat MultiDataInput SingleDataInput
% load Example1Size2.mat
%{
parfor id=1:100
    id
    [Dl,Dh]=GenerateNestedLHD(nl,nh,Dim,1e5);     
    [Dh0]=GenerateNestedLHD(nh0,nh0,Dim,1e5);     
    
    Dls(:,:,id)=Dl;
    Dhs(:,:,id)=Dh;
    Dh0s(:,:,id)=Dh0;    
end

for id =1:100
    rand
    disp(id)
    Dl=Dls(:,:,id);
    Dh=Dhs(:,:,id);
    Dh0=Dh0s(:,:,id);
    clear Yl Yh
    for jd=1:nl
        [id jd]
        Yl(jd,:)=Simulator(Dl(jd,:),1,Case);
    end
    for jd=1:nh
        Yh(jd,:)=Simulator(Dh(jd,:),2,Case);
    end
    clear  Yh0
    for jd=1:nh0
        Yh0(jd,:)=Simulator(Dh0(jd,:),2,Case);
    end
    
    MultiDataInput(id).Dl=Dl;       MultiDataInput(id).Yl=Yl;
    MultiDataInput(id).Dh= Dh;    MultiDataInput(id).Yh=Yh;
    MultiDataInput(id).XTrue=XTrue;
    MultiDataInput(id).PhysData=PhysData;    MultiDataInput(id).RatioCost=RatioCost;
    MultiDataInput(id).Budget=Budget;        MultiDataInput(id).Case=Case;
    
    SingleDataInput(id).Dl =[] ;       SingleDataInput(id).Yl=[];
    SingleDataInput(id).Dh= Dh0;    SingleDataInput(id).Yh=Yh0;
    SingleDataInput(id).XTrue=XTrue;
    SingleDataInput(id).PhysData=PhysData;    SingleDataInput(id).RatioCost=RatioCost;
    SingleDataInput(id).Budget=Budget1;        SingleDataInput(id).Case=Case;
    
end
%}
%%
%Section 2: Bayesian optimization
ZNBC_BC=1;   ZNBC_ID=0;   ZNBC_SR=2;
ZMLFSSE=1;   ZLFSSE=0; Val=1; percentage=0.99; 
for id=1:100
    id
    [T_MBC_AGP{id,1},Data_MBC_AGP{id,1}] =CalibrationAGP(MultiDataInput(id),ZNBC_BC,ZMLFSSE,Val); 'MBC-AGP'
    [T_BC_AGP{id,1},Data_BC_AGP{id,1}] =CalibrationAGP(MultiDataInput(id),ZNBC_BC,ZLFSSE,Val); 'BC-AGP'
    [T_MID_AGP{id,1},Data_MID_AGP{id,1}] =CalibrationAGP(MultiDataInput(id),ZNBC_ID,ZMLFSSE,Val); 'MID-AGP'
    [T_SR_AGP{id,1},Data_SR_AGP{id,1}] =CalibrationAGP(MultiDataInput(id),ZNBC_SR,ZLFSSE,Val); 'SR-AGP'
    [T_Nested{id,1},Data_Nested{id,1}] =CalibrationNested(MultiDataInput(id),Val); 'Nested'
    [T_SVDAGP{id,1},Data_SVDAGP{id,1}] =CalibrationSVDAGP(MultiDataInput(id),Val,percentage);'SVD-AGP'
    [T_BC_GP{id,1},Data_BC_GP{id,1}] =CalibrationBCGP(SingleDataInput(id),Val); 'BC-GP'
    [T_SR_GP{id,1},Data_SR_GP{id,1}] =CalibrationSRGP(SingleDataInput(id),Val); 'SR-GP'
    [T_SVD{id,1},Data_SVD{id,1}] =CalibrationSVD(SingleDataInput(id),Val,percentage);'SVD'
    save Example1Size2.mat
end
%%
%Section 3: Show BO results
load Example1Size2.mat
idx=(1:100);
Labels={'MBC-AGP','BC-AGP','MID-AGP','SR-AGP', 'Nested','SVD-AGP','BC-GP','SR-GP' ,'SVD'}' ;
RecordTableShow=[T_MBC_AGP(idx)  T_BC_AGP(idx)   T_MID_AGP(idx)  T_SR_AGP(idx)   T_Nested(idx)  T_SVDAGP(idx) T_BC_GP(idx)  T_SR_GP(idx)              T_SVD(idx)     ];
RecordDataShow=[Data_MBC_AGP(idx)  Data_BC_AGP(idx)  Data_MID_AGP(idx)  Data_SR_AGP(idx)    Data_Nested(idx) Data_SVDAGP(idx)   Data_BC_GP(idx)  Data_SR_GP(idx)   Data_SVD(idx)   ];

XMLE= [0.458063954892366            0.971435546875         0.999999502430791];
SSE_XMLE =3.52649277295353; 

SimMin= [ 0.12  0.0004   0];
SimMax= [  0.3   0.001   0.975];
SimGap=SimMax-SimMin;
SimXMLE=SimMin + SimGap.*XMLE;

for idxMethod=9:-1:1
    AllXhats=[];
    for idxTrain=1:numel(idx)
        Table=RecordTableShow{idxTrain,idxMethod};
        
        SSETrue_XhatsEnd(idxTrain,idxMethod)=Table.SSETrue_Xhats(end,:);
        XhatsEnd = Table.Xhats(end,:);
        L2End(idxTrain,idxMethod)=norm(XhatsEnd-XMLE);L2End_a2=L2End;
        if idxMethod<=2 || idxMethod==7
            phiEnd(idxTrain,idxMethod)=Table.phis(end,:);
        end
        
        costs=[1 RatioCost]';
        SSETrue_Xhats_iter=Table.SSETrue_Xhats;
        Xhats_iter=Table.Xhats;
        AllXhats=[AllXhats;];
        L2s_iter=sum((Xhats_iter-XMLE).^2,2).^0.5;
        Level_iter=Table.Level;
        Budget_iter=cumsum(costs(Level_iter));
        
        if idxMethod==7 || idxMethod==8 ||idxMethod==9
            TrueSSE_Xhats_Budget(1:Budget1,idxMethod,idxTrain) = interp1(Budget_iter,SSETrue_Xhats_iter,1:Budget1);
            L2s_Xhats_Budget(1:Budget1,idxMethod,idxTrain)=interp1(Budget_iter,L2s_iter,1:Budget1);
        elseif idxMethod<=4
            TrueSSE_Xhats_Budget(1:Budget,idxMethod,idxTrain) = interp1(Budget_iter,SSETrue_Xhats_iter,1:Budget);
            L2s_Xhats_Budget(1:Budget,idxMethod,idxTrain)=interp1(Budget_iter,L2s_iter,1:Budget);
        elseif idxMethod==5 || idxMethod==6
            deleteLFidx=(nl+nh+1):2:size(Table,1);
            Budget_iter(deleteLFidx,:)=[];
            SSETrue_Xhats_iter(deleteLFidx,:)=[];
            L2s_iter(deleteLFidx,:)=[];
            
            TrueSSE_Xhats_Budget(1:Budget,idxMethod,idxTrain) = interp1(Budget_iter,SSETrue_Xhats_iter,1:Budget);
            L2s_Xhats_Budget(1:Budget,idxMethod,idxTrain) = interp1(Budget_iter,L2s_iter,1:Budget);
            
        end
    end
end
meanTrueSSE_Xhats_Budget=mean(TrueSSE_Xhats_Budget(:,:,:),3);
meanTrueSSE_Xhats_Budget_SSEXMLE=meanTrueSSE_Xhats_Budget-SSE_XMLE;
meanL2s_Xhats_Budget=mean(L2s_Xhats_Budget,3);

idx1=1;
for idx2=1:9
        [ ~, ttest_p_Sh(idx2,1)]=ttest(SSETrue_XhatsEnd(:,idx1),SSETrue_XhatsEnd(:,idx2));
        [ ~, ttest_p_L2(idx2,1)]=ttest(L2End(:,idx1),L2End(:,idx2));
end
Labels1={'(i) vs (i) ','(i) vs BC-AGP ','(i) vs MID-AGP ','(i) vs SR-AGP ' ' (i) vs Nested' ,'(i) vs SVDAGP', ' (i) vs BC-GP',' (i) vs SR-GP',' (i) vs SVD'}';
Table1 =table(Labels1,mean(SSETrue_XhatsEnd)',ttest_p_Sh,mean(L2End)',ttest_p_L2)
htmlGray = [128 128 128]/255;
htmlGreen =[0.4660 0.6740 0.1880];


figure,clf
subplot(121)
FontSize6=28;
linewidth=0.1*30;
MarkerSize1=15;
plot(1:Budget,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget,1),'ko-','linewidth',linewidth+1,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget:6:Budget Budget]),hold on
plot(1:Budget,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget,2),'b:o','linewidth',linewidth+1,'MarkerSize',MarkerSize1,'MarkerFaceColor','b','MarkerIndices',[(InitialBudget+3):6:Budget Budget])
plot(1:Budget,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget,3),'k^-','linewidth',linewidth+1,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget+2:3:Budget Budget]),hold on
plot(1:Budget,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget,4),'--v','linewidth',linewidth+1,'Color', htmlGray,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget+2:3:Budget Budget])
plot(1:Budget,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget,5),':s','linewidth',linewidth+1,'color',htmlGreen,'MarkerFaceColor',htmlGreen,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget (InitialBudget):4:Budget Budget ])
plot(1:Budget,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget,6),'b:x','linewidth',linewidth+1,'MarkerSize',MarkerSize1+10,'MarkerIndices',[InitialBudget (InitialBudget+2):4:Budget Budget]),hold on
plot(1:Budget1,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget1,7),'--rs','linewidth',linewidth+1,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget InitialBudget+1:4:Budget1 Budget1])
plot(1:Budget1,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget1,8),'--h','linewidth',linewidth+1,'MarkerFaceColor','none','MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget InitialBudget+1:4:Budget1 Budget1])
plot(1:Budget1,meanTrueSSE_Xhats_Budget_SSEXMLE(1:Budget1,9),':d','linewidth',linewidth+1,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget (InitialBudget):4:Budget1 Budget1])
xlabel('Computational cost','FontWeight','normal')
ylabel('Average  $S_h(\hat{\textbf{x}}^*_{\mathbf{ML}})$-3.526493  {   }','Interpreter','latex');
leg = legend(Labels,'NumColumns',3,'Location','northeast','fontsize',20);
leg.ItemTokenSize = [64,50];
set(gca,'YScale','log','FontSize',FontSize6,'FontWeight','Bold', 'LineWidth', 3);
yticks([.0625 0.125 0.25 0.5 1 2 4 8 ])
ylim([0.05 5])
set(gca,'Position',[0.1 0.2 0.41 0.71])
title('(a)','FontSize',25,'FontWeight','Bold')
xticks([InitialBudget:2:Budget1 ])
xlim([InitialBudget-0.1 Budget1+0.2])

subplot(122)
plot(1:Budget,meanL2s_Xhats_Budget(1:Budget,1),'ko-','linewidth',linewidth+1,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget:6:Budget Budget]),hold on
plot(1:Budget,meanL2s_Xhats_Budget(1:Budget,2),'b:o','linewidth',linewidth+1,'MarkerSize',MarkerSize1,'MarkerFaceColor','b','MarkerIndices',[(InitialBudget+3):6:Budget Budget])
plot(1:Budget,meanL2s_Xhats_Budget(1:Budget,3),'k^-','linewidth',linewidth+1,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget+2:3:Budget Budget]),hold on
plot(1:Budget,meanL2s_Xhats_Budget(1:Budget,4),'--v','linewidth',linewidth+1,'Color', htmlGray,'MarkerSize',MarkerSize1,'MarkerIndices',[(InitialBudget+2):3:Budget Budget])
plot(1:Budget,meanL2s_Xhats_Budget(1:Budget,5),'g:s','linewidth',linewidth+1,'color',htmlGreen,'MarkerFaceColor',htmlGreen,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget (InitialBudget):4:Budget Budget ])
plot(1:Budget,meanL2s_Xhats_Budget(1:Budget,6),'b:x','linewidth',linewidth+1,'MarkerSize',MarkerSize1+10,'MarkerIndices',[InitialBudget (InitialBudget):3:Budget Budget]),hold on
plot(1:Budget1,meanL2s_Xhats_Budget(1:Budget1,7),'--s','linewidth',linewidth+1,'Color', 'r','MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget InitialBudget+1:4:Budget1 Budget1])
plot(1:Budget1,meanL2s_Xhats_Budget(1:Budget1,8),'--h','linewidth',linewidth+1,'MarkerFaceColor','none','MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget InitialBudget+1:4:Budget1 Budget1])
plot(1:Budget1,meanL2s_Xhats_Budget(1:Budget1,9),':d','linewidth',linewidth+1,'MarkerSize',MarkerSize1,'MarkerIndices',[InitialBudget (InitialBudget):4:Budget1 Budget1 ])
xlim([InitialBudget,Budget+1])
xlabel('Computational cost','FontWeight','normal')
ylabel('Average  $L_2(\hat{\textbf{x}}^*_{\mathbf{ML}})$','Interpreter','latex');  
yticks([0.2 0.3:0.1:1.2 ] )
leg = legend(Labels,'NumColumns',3,'Location','northeast','fontsize',20);
leg.ItemTokenSize = [64,45];
set(gca,'Position',[0.583 0.2 0.41 0.71])
title('(b)','FontSize',24,'FontWeight','Bold')
ylim([0.19  0.83])
set(findobj(gcf,'type','axes'),'FontSize',FontSize6,'FontWeight','Bold', 'LineWidth', 3);
set(gcf,'Position',[          0         100        1920         615])
xticks([InitialBudget:2:Budget1 ])
xlim([InitialBudget-0.1 Budget1+0.2])


































