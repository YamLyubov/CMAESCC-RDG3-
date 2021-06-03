
function group = grouping(seps, nonseps)
    
    % further divide the separable subcomponents
    div_sepsize = 100;
    div_num = ceil(length(seps)/div_sepsize);
    div_seps = cell(1, div_num);
    div_startpts = 1:div_sepsize:length(seps);
    for i = 1:div_num-1
       div_seps{i} = seps(div_startpts(i):div_startpts(i+1)-1);
    end
    if (~isempty(seps))
       div_seps{div_num} = seps(div_startpts(div_num):end);
    end
    group = [div_seps nonseps];
end
