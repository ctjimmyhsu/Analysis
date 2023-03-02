%% THIS PROGRAM IS USED TO ANALYZE THE MT SIMULATION FOR GROWTH RATE, LIFETIME AND SHRINKAGE RATE
function [mgrate,mlifetime,mshkrate,mpause]...
                    = MTParaSearch_anal_new(data,tseries,closedring,sim_num,uppercut,lowercut)
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
alpha = 0.1;       %Percentage of slope difference
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
    idxmin = 0;
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
                if abs(minpks(j,1))<lowercut && maxpks(j,1)>uppercut
                    %% GROWTH PART
                    idxmax = find(Ldata(:,1)==maxlocs(j,1));
                    idxmin = find(Ldata(:,1)==minlocs(j,1));
                    dmax = ceil(abs(idxmax-idxmin).*0.1);
                    dmin = dmax;
                    if (idxmin-dmin)<1
                        idmi = 1;
                    else
                        idmi = idxmin-dmin;
                    end
                    [M,I] = min(Ldata(idmi:idxmin,2));
                    newmin = idmi+I-1;%Shift to correct local min
                    %minlocs(j,1) = Ldata(newmin,1); %update value
                    %minpks(j,1) = -Ldata(newmin,2); %update value
                    if (idxmax+dmax)>l(1,1)
                        idmx = l(1,1);
                    else
                        idmx = (idxmax+dmax);
                    end
                    [M,I] = max(Ldata(idxmax:idmx,2));
                    newmax = idxmax+I-1;%shift to correct local maximum
                    %maxlocs(j,1)= Ldata(newmax,1); %update value
                    %maxpks(j,1) = Ldata(newmax,2); %update value
                    div = abs(Ldata(newmax,2)-Ldata(newmin,2)).*alpha./2;
                    %div1 = abs((Ldata(newmax,2)-Ldata(newmin,2))./(newmax-newmin));
                    P = polyfit(Ldata(newmin:newmax,1),Ldata(newmin:newmax,2),1);
                    [p1,p2] = RANSAC(Ldata(newmin:newmax,1:2)',2,1000,div,0.5);
                    if p1>0
                        grate(gcount,1) = p1.*60;
                        gcount = gcount+1;
                        t1 = Ldata(newmin,1);
                        t2 = Ldata(newmax,1);
                        lifetime(lcount,1) = t2-t1;
                        lcount = lcount+1;
                    end  
                    %% SHRINKAGE PART
                    if ((j+1)<=minlength(1,1) && abs(minpks(j+1,1))<lowercut)%check next min for shrinkage
                        nangle = angle(:,1).*(abs(angle(:,1))>0.079829985);%Angle that is bigger than adding in 1 dimer
                        idxmax = newmax;%The new maximum from above
                        idxmin = find(Ldata(:,1)==minlocs(j+1,1));
                        dmin = ceil(abs(idxmin-idxmax).*0.4);
                        if (idxmin-dmin)<1
                            idmi = 1;
                        else
                            idmi = idxmin-dmin;
                        end
                        [M,I] = min(Ldata(idmi:idxmin,2));
                        newmin2 = idmi+I-1;%Shift to new local minimum
                        %minlocs(j+1,1) = Ldata(newmin2,1); %update value
                        %minpks(j+1,1) = -Ldata(newmin2,2); %update value
                        [Ma,Ia] = min(nangle(idxmax:newmin2,1));%Find the maximum slope of <L>
                        maxangle = idxmax+Ia-1;
                        [Mc,Ic] = min(cangle(idxmax:newmin2,1));%Find the maximum slope closedring
                        closedangle = idxmax+Ic-1;
                        for k=maxangle:-1:idxmax %Goping towards the left
                            if (nangle(k-2,1)==0&&nangle(k-1,1)==0&&nangle(k,1)==0)
                                if (k<closedangle)%left of the closed angle
                                    shkstart = k;
                                    break;
                                end
                            end
                        end
                        pause(pcount,1) = abs(Ldata(idxmax,1)-Ldata(shkstart,1));
                        pcount = pcount+1;
                        %maxangle
                        %newmin2
                        for k=maxangle:newmin2 %Goping towards the right
                            if (nangle(k,1)==0&&nangle(k+1,1)==0&&nangle(k+2,1)==0)
                                if (k>closedangle&&Ldata(k,2)<lowercut)%left of the closed angle
                                    shkend = k;
                                    break;
                                end
                            end
                        end
                        div = abs(Ldata(shkstart,2)-Ldata(shkend,2)).*alpha./2;
                        %div1 = abs((Ldata(shkstart,2)-Ldata(shkend,2))./(shkend-shkstart));
                        %shkstart
                        %shkend
                        P = polyfit(Ldata(shkstart:shkend,1),Ldata(shkstart:shkend,2),1);
                        [p1,p2] = RANSAC(Ldata(shkstart:shkend,1:2)',2,100,div,0.5);
                        if(p1<0)%Only if it's a negative slope
                            shkrate(scount,1) = abs(p1*60);
                            scount = scount+1;
                        end
                    end%End shrinkage part
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