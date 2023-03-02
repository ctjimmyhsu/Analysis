%% THIS PROGRAM IS USED TO ANALYZE THE MT SIMULATION FOR GROWTH RATE, LIFETIME AND SHRINKAGE RATE
function [mgrate,mlifetime,mshkrate,mpause]...
                    = MTParaSearch_anal_new2(data,tseries,closedring,sim_num,uppercut,lowercut)
%% INITIALIZATION
gcount = 1;%Counts for grate
lcount = 1;%Counts for lifetime
scount = 1;%Counts for shkrate
pcount = 1;%Counts for pauses
grate = zeros(5,1);%grate
lifetime = zeros(5,1);%lifetime
shkrate = zeros(5,1);%shkrate
pause = zeros(5,1);%pauses
l = size(data);     %size of full matrix
alpha = 0.1;        %Percentage of slope difference
dt=0.1;             %Time steps for MT trace
for i=1:sim_num
    %% SMOOTH THE DATAE
    time = tseries(1:15:end,i);
    yy2 = smooth(time,data(1:15:end,i),'lowess');
    yy2(1,1) = 0;%Reset to 0 incase of smoothing
    [maxpks,maxlocs] = findpeaks(yy2,time);  %Find local max
    [minpks,minlocs] = findpeaks(-yy2,time); %Find local min
    Ldata = horzcat(tseries(:,i),data(:,i)); %For RANSAC 2D array
    Cdata = horzcat(tseries(:,i),closedring(:,i));%For correcting
    angle = atan((Ldata(2:end,2)-Ldata(1:end-1,2))./(Ldata(2:end,1)-Ldata(1:end-1,1)));%Angle between points
    cangle = atan((Cdata(2:end,2)-Cdata(1:end-1,2))./(Cdata(2:end,1)-Cdata(1:end-1,1)));
    idxmax = 0;
    idxmin1 = 0;
    %Note pks = N by 1 array, locs = N by 1 array
    %% CALCULATION
    if isempty(maxlocs)~=1 && isempty(minlocs)~=1 %the MT fall apart in this trace
        minpks = vertcat(0,minpks);     %artificial minimum of origin
        minlocs = vertcat(0,minlocs);   %artificial minimum of origin
        maxlength = size(maxlocs);
        minlength = size(minlocs);
        %% LOCAL MINIMUM FIRST
        if minlocs(1,1) < maxlocs(1,1)
            if maxlength(1,1)>minlength(1,1)
                length = minlength(1,1);
            elseif maxlength(1,1)<minlength(1,1)
                length = maxlength(1,1);
            elseif maxlength(1,1)==minlength(1,1)
                length = minlength(1,1);
            end
            for j=1:length %%Go through all peak and trough
                if abs(minpks(j,1))<lowercut && maxpks(j,1)>uppercut && (j+1)<=minlength(1,1) && abs(minpks(j+1,1))<lowercut
                    %% GETTING THE CORRECT INDICES
                    idxmin1 = find(Ldata(:,1)==minlocs(j,1));   %First min
                    idxmax = find(Ldata(:,1)==maxlocs(j,1));    %Max
                    idxmin2 = find(Ldata(:,1)==minlocs(j+1,1)); %Second min
                    %dg = ceil(abs(idxmax-idxmin1).*0.1);    %Points search for max
                    %ds = ceil(abs(idxmin2-idxmax).*0.4);    %Points search for min
                    dg = 4/dt;                              %Points search +/- 4s
                    ds = 4/dt;
                    if (idxmin1==1)%First point
                        newmin1 = idxmin1;
                    else
                        [M,I] = min(Ldata((idxmin1-ds):idxmin1,2));         %Search for true min
                        newmin1 = (idxmin1-ds)+I-1;  %Shift to new local minimum
                    end
                    [M,I] = max(Ldata((idxmax-dg):(idxmax+dg),2));  %Search for true max
                    newmax = (idxmax-dg)+I-1;   %Shift to new local maxmimum
                    [M,I] = min(Ldata((idxmin2-ds):idxmin2,2));         %Search for true min
                    newmin2 = (idxmin2-ds)+I-1;  %Shift to new local minimum
                    %% FITTING SECTION: GROWTH
                    divg = abs(Ldata(newmax,2)-Ldata(newmin1,2)).*alpha./2;
                    [g1,g2] = RANSAC(Ldata(newmin1:newmax,1:2)',2,1000,divg,0.5);
                    if g1>0
                        gtemp = g1; %For self contained
                        grate(gcount,1) = g1.*60;
                        gcount = gcount+1;
                        t1 = Ldata(newmin1,1);
                        t2 = Ldata(newmax,1);
                        lifetime(lcount,1) = t2-t1;
                        lcount = lcount+1;
                    else
                        G = polyfit(Ldata(newmin1:newmax,1),Ldata(newmin1:newmax,2),1);
                        gtemp = G(1,1);
                        grate(gcount,1) = G(1,1).*60;
                        gcount = gcount+1;
                        t1 = Ldata(newmin1,1);
                        t2 = Ldata(newmax,1);
                        lifetime(lcount,1) = t2-t1;
                        lcount = lcount+1;
                    end  
                    %% FITTING SECTION: SHRINKAGE
                    divs = abs(Ldata(newmax,2)-Ldata(newmin2,2)).*alpha./2;
                    nangle = angle(:,1).*(abs(angle(:,1))>0.079829985);%Angle not bigger than grate are all 0
                    [Ma,Ia] = min(nangle(newmax:newmin2,1));%Find the maximum slope of <L>
                    maxangle = newmax+Ia-1;
                    %[Mc,Ic] = min(cangle(newmax:newmin2,1));%Find the maximum slope closedring
                    %closedangle = newmax+Ic-1;
                    %nangle = angle(:,1).*(abs(angle(:,1))>0.079829985);%Angle that is bigger than adding in 1 dimer
                    for k=newmax:maxangle
                        if (nangle(k+1,1)<0&&nangle(k+2,1)<0&&nangle(k+3,1)<0)
                            %if (k<closedangle)%left of the closed angle
                            shkstart = k;
                            break;
                            %end
                        end
                    end
                    pause(pcount,1) = abs(Ldata(newmax,1)-Ldata(shkstart,1));
                    pcount = pcount+1;
                    %maxangle
                    %newmin2
                    for k=maxangle:newmin2 %Goping towards the right
                        if (nangle(k,1)==0&&nangle(k+1,1)==0&&nangle(k+2,1)==0)
                            %if (k>closedangle&&Ldata(k,2)<lowercut)%left of the closed angle
                            shkend = k;
                            break;
                        else
                            shkend = k;
                            %end
                        end
                    end
                    %P = polyfit(Ldata(shkstart:shkend,1),Ldata(shkstart:shkend,2),1);
                    [p1,p2] = RANSAC(Ldata(shkstart:shkend,1:2)',2,1000,divs,0.5);
                    if(p1<0)%Only if it's a negative slope
                        shkrate(scount,1) = abs(p1*60);
                        scount = scount+1;
                    else
                        S = polyfit(Ldata(shkstart:shkend,1),Ldata(shkstart:shkend,2),1);
                        shkrate(scount,1) = abs(S(1,1)*60);
                        scount = scount+1;
                    end
                    %                     figure,
                    %                     plot(Ldata(newmin:newmax,1),Ldata(newmin:newmax,2),'-rs','LineWidth',2,'MarkerFaceColor','r','MarkerSize',3);
                    %                     hold on
                    %                     plot(Ldata(newmax:shkstart,1),Ldata(newmax:shkstart,2),'-bo','LineWidth',2,'MarkerFaceColor','b','MarkerSize',3);
                    %                     plot(Ldata(shkstart:shkend,1),Ldata(shkstart:shkend,2),'-gd','LineWidth',2,'MarkerFaceColor','g','MarkerSize',3);
                    %                     hold off
                    %                     %pbaspect([1 1 1])
                    %                     xlabel('Time (s)','fontSize',20)
                    %                     ylabel('<L>','fontSize',20)
                    %                     set(gca,'fontSize',20)
                    %                     %legend('Growth','Slow shrinkage','fast Shrinkage','location','NorthWest')
                    %                     figureHandle = gcf;
                    %                     set(findall(figureHandle,'type','text'),'fontSize',20)
                end%End of qualifying growth               
            end%End of Length loop
        end
    end
end
mgrate = grate;
mlifetime = lifetime;
mshkrate = shkrate;
mpause = pause;