%%This is a code to simulate MT length using gamma distribution catastrophe
%%frequency: Recode from Sami 
%% PREPARE FOR SIMULATION
close all
clear all
Gamma = importdata('Bovine_Gamma.txt'); %Read in gamma parameter data
Growth = importdata('Bovine_grate.txt');%Read in growth rate data
time = 1200.0;    %total time of simulation [s]
dt = 5.0;         %Sampling time interval [s]
vs = 15.0;        %MT shrinkage rate [um/min]
nuclflag = 0;  %Nucleation flag 1: on, 0: off
nuctime = 50.0;  %Nucleation half time [s]
N = 1E4;        %Total number of microtubule
Lwidth = 0.5;     %Bin width
Lmax = 5;      %Max value of the bins
Lx = Lwidth/2:Lwidth:Lmax-(Lwidth/2);%x-axis for plotting
plotdt = 200;           %Output plotting 
CC = jet(time/plotdt);  %Color for plotting
count=1;

%% PARAMETER SETUP
vg = Growth(1,2);   %Growth Rate
steps = Gamma(1,2); %step parameter 
rate = Gamma(1,3)/60;  %Rate parameter: scale = 1/rate, change to seconds
shift = Gamma(1,4); %Shift from 0
lengths = zeros(N,1);   %Number of MT
states = zeros(N,1);    %State of the MT: -1=shrinking,0=nucleating,1=growing
catadeci = zeros(N,1);   %Catastrophe decisions
nucldeci = zeros(N,1);   %Nucleation decisions
nucltime = zeros(N,1);  %Tracking the time interval MT in nucleation phase
growtime = zeros(N,1);  %Tracking the time interval MT in growth phase


%% START SIMULATION
for t = 0:dt:time       %Go through simulation time
    avglength = 0;      %Average length of microtubule
    for i=1:N           %Go through number of MT
        if states(i,1) == 0   %IN START PHASE
            if nuclflag==1
                nucldeci(i,1) = exprnd(nuctime);
            else
                nucldeci(i,1) = 0;
            end
            states(i,1) = 1;%Go to growing phase
        elseif states(i,1) == 1%IN NUCLEATION PHASE
            nucltime(i,1) = nucltime(i,1) + 1;   %Add count
            currnucltime = nucltime(i,1)*dt;%Calculate time
            if currnucltime > nucldeci(i,1)
                states(i,1) = 2;    %Go to growing phase
                nucltime(i,1) = 0;
            end
        elseif states(i,1) == 2 %IN GROWTH PHASE   
            catadeci(i,1) = gamrnd(steps,(1./rate));%Find time that catastrophe
            states(i,1) = 3;        %Go to after catastrophe
        elseif states(i,1) == 3 %FIND TIME TO CATASTROPHE    
            lengths(i,1) = lengths(i,1) + (vg./60).*dt;%Unit [um]
            growtime(i,1) = growtime(i,1) + 1;
            currlifetime = growtime(i,1)*dt;
            if currlifetime > catadeci(i,1) %Passed catastrophe time
                states(i,1) = -1;   %Go to shrinkage phase
                growtime(i,1) = 0;  %Reset growth time interval
            end
        elseif states(i,1) == -1    %IN SHRINKAGE PHASE
            lengths(i,1) = lengths(i,1) - (vs./60).*dt;%Unit [um]
            if lengths(i,1) < 0
                lengths(i,1) = 0;
                states(i,1) = 0;
            end
        end%End state check        
    end%End MT loop
    if mod(t,plotdt)==0&&t~=0
        Lhist = histogram(lengths(:,1),'BinWidth',Lwidth,'BinLimits',[0 Lmax],'Normalization','probability'...
            ,'Visible','off');
        if count==1
            Lprob = Lhist.Values;
        else
            Lprob = vertcat(Lprob,Lhist.Values);
        end
%         if count ==1
%             figure, 
%             plot(Lx,Lprob,'color',CC(count,:),'LineWidth',2,'Marker','s','MarkerFaceColor',CC(count,:),'MarkerSize',3);
%             hold on
%         else
%             plot(Lx,Lprob,'color',CC(count,:),'LineWidth',2,'Marker','s','MarkerFaceColor',CC(count,:),'MarkerSize',3);
%         end
         count = count+1;
    end
end%End time loop
figure,
hold on
for i=1:(time/plotdt)
    plot(Lx,Lprob(i,:),'color',CC(i,:),'LineWidth',2,'Marker','s','MarkerFaceColor',CC(i,:),'MarkerSize',3);
end
avglength = mean(lengths(:,1));
xavg = 0:0.1:5;
p = exp(-Lx./avglength)./avglength;
z = sum(p);
plot(Lx,p./z,'-ko');
xlabel('Length (\mum)','fontSize',20)
ylabel('P(Length)','fontSize',20)
legend('t=200s','t=400s','t=600s','t=800s','t=1000s','t=1200s','exp','location','NorthEast')
set(gca,'fontSize',20)
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',20)
box on
hold off


