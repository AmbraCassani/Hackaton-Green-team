function namelist = create_list(direc_name)

direc = dir(direc_name);


namelist = cell(length(direc),1);
for i=1:length(direc)
    namelist{i,1} = direc(i).name(1:end-4); % 1:10 gives only the name without the format, usefull for the main for loop
end

end