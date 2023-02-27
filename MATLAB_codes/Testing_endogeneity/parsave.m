function parsave(fname, x)
%     if exist(fname, 'file')==2
%       %name = split(fname,'.mat');
%       save(fname, 'x', '-append')
%     else
        save(fname, 'x')
%     end
end
