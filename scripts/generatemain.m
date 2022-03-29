% Read txt into cell A
fid = fopen('test_tickness.m','r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);
ints=1;
mkdir(['main' num2str(ints)])
for tt=1:12
    name=(tt-1)*(10^2)*7;
    for mm=1:1:7
        for kd=1:1:5
            for kde=1:1:5
                
                name=name+1;
                B=A;
                % Change cell A
                B{7} =  strrep(A{7},'(1)',['(',num2str(tt),')']);
                B{8} =  strrep(A{8},'(1)',['(',num2str(mm),')']);
                B{9} =  strrep(A{9},'(1)',['(',num2str(kd),')']);
                B{10} =  strrep(A{10},'(1)',['(',num2str(kde),')']);
                
                % Write cell A into txt
                fid = fopen(['main' num2str(ints) '/main' num2str(name) '.m'], 'w');
                for i = 1:numel(B)
                    if B{i+1} == -1
                        fprintf(fid,'%s', B{i});
                        break
                    else
                        fprintf(fid,'%s\n', B{i});
                    end
                end
                fclose(fid);
            end
        end
    end
end

