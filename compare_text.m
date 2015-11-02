function similarity=compare_text(n,m,dirx)
%
%--------------------------------------------------------------------------
% Need to set the directory
cd(dirx); 
struct_dir = dir('*.txt');
for i = 1:length(struct_dir)
    if struct_dir(i).isdir == 0
        fidx = fopen(struct_dir(i).name,'r');
        textcell = textscan(fidx,'%s')
        fclose(fidx);
    end
end