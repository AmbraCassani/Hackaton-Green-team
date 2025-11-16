function [struct1,struct2,idx_co,idx_fl] = memorization_subgroup(struct1,struct2,namelist,total_walking_percent, total_steps, avgstride,i,idx_co,idx_fl)

% Similarly to the other "memorization function" it stores collected data
% into final structures



if strcmp(namelist{i}(1:2), 'CO')


    struct1.walking_percent(idx_co) = total_walking_percent;
    struct1.total_steps(idx_co) = total_steps;
    struct1.avgstride(idx_co) = avgstride;



    idx_co = idx_co +1;

else
    struct2.walking_percent(idx_fl) = total_walking_percent;
    struct2.total_steps(idx_fl) = total_steps;
    struct2.avgstride(idx_fl) = avgstride;



    idx_fl = idx_fl +1;
end


