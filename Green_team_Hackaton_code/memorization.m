function [struct1,struct2,idx_co,idx_fl] = memorization(struct1,struct2,namelist,alpha,h_zejl,t_zejl,i,idx_co,idx_fl,sgl,idx_turns,fs,pos_peak,win_sma)

% The function takes as input the structures in which the main results of
% the gait speed analysis will be stored, together with other parameters
% helpfull for the following of the analysis.


if strcmp(namelist{i}(1:2), 'co')
    struct1.nsubject = idx_co;
    mean_height = 1.64; %meters, given by the main paper
    leg_length = alpha*mean_height; % through the function parameters alpha it is possible to change this relationship
    step_length(i) = 2*sqrt(2*leg_length*h_zejl(i) - h_zejl(i)^2); % inverted pendulum model - Zejlstra
    struct1.gaitspeedZejlstra(idx_co) = step_length(i) * length(pos_peak)/t_zejl(end); 




    BA = find_BA(sgl,idx_turns,fs,1);
    SMAc = movmean( abs(BA(:,1))+ abs(BA(:,2)) + abs(BA(:,3)), win_sma);
    struct1.SMA(idx_co) = mean(SMAc); % the SMA is in this way stored and will be used after as a threshold



    idx_co = idx_co +1;

else
    struct2.nsubject = idx_fl;
    mean_height = 1.61; %meters, given by the main paper
    leg_length = alpha*mean_height; 
    step_length(i) = 2*sqrt(2*leg_length*h_zejl(i) - h_zejl(i)^2);
    struct2.gaitspeedZejlstra(idx_fl) = step_length(i) * length(pos_peak)/t_zejl(end);



    BA = find_BA(sgl,idx_turns,fs,1);
    SMAf = movmean( abs(BA(:,1))+ abs(BA(:,2)) + abs(BA(:,3)), win_sma);
    struct2.SMA(idx_fl) = mean(SMAf);


    idx_fl = idx_fl +1;



end
