%% THIS PROGRAM IS USED TO ANALYZE THE MT SIMULATION FOR GROWTH RATE, LIFETIME AND SHRINKAGE RATE
function [mgrate,stdgrate,maxgrate,mingrate,lgrate,mlifetime,stdlifetime,maxlifetime,minlifetime,...
    llifetime,mshkrate,stdshkrate,maxshkrate,minshkrate,lshkrate]...
                    = MTParaSearch_anal(data,tseries,closedring,sim_num,cutoff)
%% INITIALIZATION
gcount = 1;%Counts for grate
lcount = 1;%Counts for lifetime
scount = 1;%Counts for shkrate
grate = zeros(5,1);%grate
lifetime = zeros(5,1);%lifetime
shkrate = zeros(5,1);%shkrate
l = size(data);     %size of full matrix
dd=30;
for i=1:sim_num
    %% SMOOTH THE DATAE
    time = tseries(1:15:end,i);
    yy2 = smooth(time,data(1:15:end,i),'lowess');
    yy2(1,1) = 0;%Reset to 0 incase of smoothing
    [maxpks,maxlocs] = findpeaks(yy2,time);  %Find local max
    [minpks,minlocs] = findpeaks(-yy2,time); %Find local min
    tdata = horzcat(tseries(:,i),data(:,i)); %For RANSAC 2D array
    adata = horzcat(tseries(:,i),closedring(:,i));
    angle = atan((adata(2:end,2)-adata(1:end-1,2))./(adata(2:end,1)-adata(1:end-1,1)));
    t1 = 0;
    t2 = 0;
    %Note pks = N by 1 array, locs = N by 1 array
    %% CALCULATION    
    if isempty(maxlocs)~=1 && isempty(minlocs)~=1
        maxlength = size(maxlocs);
        minlength = size(minlocs);      
        if minlocs(1,1) < maxlocs(1,1)
            %% LOCAL MINIMUM FIRST
            if maxlength(1,1)>minlength(1,1)
                length = minlength(1,1);
            elseif maxlength(1,1)<minlength(1,1)
                length = maxlength(1,1);
            elseif maxlength(1,1)==minlength(1,1)
                length = minlength(1,1);
            end
            for j=1:length
                if abs(minpks(j,1))<cutoff && maxpks(j,1)>cutoff
                    %Growth Rate
                    idxmax = find(tdata(:,1)==maxlocs(j,1));
                    idxmin = find(tdata(:,1)==minlocs(j,1));
                    %di = abs(idxmax-idxmin);
                    %dd = round(di*0.05);
                    if (idxmin-dd)<1
                        idmi = 1;
                    else
                        idmi = idxmin-dd;
                    end
                    [M,I] = min(tdata(idmi:idxmin,2));
                    In = idmi+I-1;%Shift to correct local min
                    if (idxmax+dd)>l(1,1)
                        idmx = l(1,1);
                    else
                        idmx = (idxmax+dd);
                    end
                    [M,I] = max(tdata(idxmax:idmx,2));
                    Ix = idxmax+I-1;%shift to correct local maximum
                    [p1,p2] = RANSAC(tdata(In:Ix,1:2)',2,1000,0.008,0.5);
                    if p1>0
                        grate(gcount,1) = p1.*60;
                        gcount = gcount+1;
                        t1 = tdata(In,1);
                        t2 = tdata(Ix,1);
                        lifetime(lcount,1) = t2-t1;
                        lcount = lcount+1;
                    end
                end
            end
            %Shrinkage rate
            if minlength(1,1)>1%check first value
                minpks = minpks(2:end,1);
                minlocs = minlocs(2:end,1);
                length = size(minlocs);
                for j=1:length(1,1)
                    if maxpks(j,1)>cutoff && abs(minpks(j,1))<cutoff  
                        idxmax = find(tdata(:,1)==maxlocs(j,1));
                        idxmin = find(tdata(:,1)==minlocs(j,1));
                        %di = abs(idxmax-idxmin);
                        %dd = round(di*0.05);%Find Max
%                         if (idxmax+dd)>l(1,1)
%                             idmx = l(1,1);
%                         else
%                             idmx = (idxmax+dd);
%                         end
%                         [M,I] = max(tdata(idxmax:idmx,2));
%                         Ix = idxmax+I-1;%shift to correct local maximum
                        if (idxmin-dd)<1
                            idmi = 1;
                        else
                            idmi = idxmin-dd;
                        end
                        [M,I] = min(tdata(idmi:idxmin,2));
                        In = idmi+I-1;%Shift to correct local min
                        [Ma,Ia] = min(angle(idxmax:idxmin));%Find the maximum slope
                        L = In-(idxmax+Ia-1);%Distance between the min and max slope
                        [p1,p2] = RANSAC(tdata(In-(2*L):In,1:2)',2,100,0.05,0.5);
                        %[p1,p2] = RANSAC(tdata(idxmax:idxmin,1:2)',2,100,0.008,0.1);
                        shkrate(scount,1) = p1*60;
                        scount = scount+1;
                        %shkrate(scount,1) = abs((maxpks(j,1)-abs(minpks(j,1)))./(maxlocs(j,1)-minlocs(j,1)).*60);
                        %scount = scount+1;
                    end
                end
            end           
        elseif minlocs(1,1) > maxlocs(1,1)&&minlength(1,1)>1 
            %% LOCAL MAXMIMUM FIRST
            %Shrinkage rate
            if maxlength(1,1)>minlength(1,1)
                length = minlength(1,1);
            elseif maxlength(1,1)<minlength(1,1)
                length = maxlength(1,1);
            elseif maxlength(1,1)==minlength(1,1)
                length = minlength(1,1);
            end
            for j=1:length
                if maxpks(j,1)>cutoff && abs(minpks(j,1))<cutoff
                    idxmax = find(tdata(:,1)==maxlocs(j,1));
                    idxmin = find(tdata(:,1)==minlocs(j,1));
                    di = abs(idxmax-idxmin);
                    ddd = round(di*0.4);
                    %                     if (idxmax+dd)>l(1,1)%Find Max
                    %                         idmx = l(1,1);
                    %                     else
                    %                         idmx = (idxmax+dd);
                    %                     end
                    %                     [M,I] = max(tdata(idxmax:idmx,2));
                    %                    Ix = idxmax+I-1;%shift to correct local maximum
                    if (idxmin-ddd)<1%Find Min
                        idmi = 1;
                    else
                        idmi = idxmin-ddd;
                    end
                    [M,I] = min(tdata(idmi:idxmin,2));
                    In = idmi+I-1;%Shift to correct local min
                    [Ma,Ia] = min(angle(idxmax:idxmin));%Find the maximum slope
                    L = In-(idxmax+Ia-1);%Distance between the min and max slope
                    %P = polyfit(tdata(In-(L):In,1),tdata(In-(L):In,2),1);
                    if L>0
                        [p1,p2] = RANSAC(tdata(In-(L):In,1:2)',2,100,0.05,0.99);
                        %p1 = P(1,1);
                        if(p1<0)%Only if it's a negative slope
                            shkrate(scount,1) = abs(p1.*60);
                            scount = scount+1;
                        end
                    end
                end
            end
            maxpks = maxpks(2:end);
            maxlocs = maxlocs(2:end);
            length = size(maxlocs);
            for j=1:length(1,1)
                if abs(minpks(j,1))<cutoff && maxpks(j,1)>cutoff
                    %Grate 
                    idxmin = find(tdata(:,1)==minlocs(j,1));
                    idxmax = find(tdata(:,1)==maxlocs(j,1));
                    di = abs(idxmax-idxmin);
                    %dd = round(di*0.05);
                    if (idxmin-dd)<1
                        idmi = 1;
                    else
                        idmi = idxmin-dd;
                    end
                    [M,I] = min(tdata(idmi:idxmin,2));
                    In = idmi+I-1;%Shift to correct local min
                    if (idxmax+dd)>l(1,1)
                        idmx = l(1,1);
                    else
                        idmx = (idxmax+dd);
                    end
                    [M,I] = max(tdata(idxmax:idmx,2));
                    Ix = idxmax+I-1;%shift to correct local maximum
                    [p1,p2] = RANSAC(tdata(In:Ix,1:2)',2,100,0.008,0.3);
                    if p1>0
                        grate(gcount,1) = p1.*60;
                        gcount = gcount+1;
                        t1 = tdata(In,1);
                        t2 = tdata(Ix,1);
                        lifetime(lcount,1) = t2-t1;
                        lcount = lcount+1;
                    end
                end
            end
        end
%     elseif isempty(maxlocs)==1 && isempty(minlocs)==1 %only growing
%         grate(gcount,1) = (yy2(end-5,1)-yy2(5,1))./(time(end-5,1)-time(5,1)).*60;
%         gcount = gcount+1;
     end
end
mgrate = mean(grate);
stdgrate = std(grate);
maxgrate = max(grate);
mingrate = min(grate);
lgrate = max(size(grate));
mlifetime = mean(lifetime);
stdlifetime = std(lifetime);
maxlifetime = max(lifetime);
minlifetime = min(lifetime);
llifetime = max(size(lifetime));
mshkrate = shkrate;
stdshkrate = std(shkrate);
maxshkrate = max(shkrate);
minshkrate = min(shkrate);
lshkrate = max(size(shkrate));
% mgrate = grate;
% stdgrate = std(grate);
% mlifetime = lifetime;
% stdlifetime = std(lifetime);
% mshkrate = shkrate;
% stdshkrate = std(shkrate);