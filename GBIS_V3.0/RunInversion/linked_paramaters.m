function [m] = linked_paramaters(m,link)
% Function to check if any parameters need to be held the same

n_param=size(m,1);
link_param=find(link);

if ~isempty(link_param)
    for ii=1:numel(link_param)
        m(link_param(ii))=m(link(link_param(ii)));
    end
end

