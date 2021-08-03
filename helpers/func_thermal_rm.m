function [t_out,data_out,fit_out] = func_thermal_rm(t_in,data_in,order,startval,endval) 
t2=t_in(startval:endval);

data_out = zeros(size(data_in,1),size(data_in,2),length(startval:endval));

for x = 1:size(data_in,1);
    for y = 1:size(data_in,2);
        data_tmp = squeeze(data_in(x,y,startval:endval));
        [p,s,mu] = polyfit(t2',data_tmp,order);
        fit = polyval(p,t2',[],mu);
        data_out(x,y,:)= data_tmp-fit;
        
        if nargout == 3
            fit_out(x,y,:) = fit;
        end
    end
end
t_out=t2;