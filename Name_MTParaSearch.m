%% This function creates a name for data input for parameter search
function [ datastr ] = Name_MTParaSearch(str1,prn,paraval,mfactor,con,tpoint)
for i=1:prn
    num =num2str(abs(paraval(1,i).*mfactor(1,i)));
    if i==1;
        datastr = strcat(str1,num,'_');
    elseif i==prn
        datastr = strcat(datastr,num);
    else
        datastr = strcat(datastr,num,'_');
    end
end
datastr = strcat(datastr,'_',num2str(con));
datastr = strcat(datastr,'_',num2str(tpoint));