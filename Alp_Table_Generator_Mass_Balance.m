BEM = 9165;
BEM_arm = 292.18;
datasheet = 'Post_Flight_Datasheet_Flight_2_DD_14_3_2018.xlsx';

M = Mass_CG_Moment(BEM, BEM_arm, datasheet);    %[kg, m, Nm]

%Set number of digits that you want here!!!!!!!
digits_Mass = 2;
digits_CG = 4;
digits_Moment = 1;

%Don't look too much into what textarray does, I pretty much forgot myself
%too
fid = fopen('Mass_Balance_Table.txt','w');
textarray = string(num2str(M(:,1),char('%.'+string(digits_Mass)+'f')))+' & '+...
    string(num2str(M(:,2),char('%.'+string(digits_CG)+'f')))+' & '+...
    string(num2str(M(:,3),char('%.'+string(digits_Moment)+'f')));

text = textarray(1);
for i = 2:length(textarray)
    text = text + '\r\n' + textarray(i);
end

%The table file is written and saved
fprintf(fid,text);
fclose(fid);

!Mass_Balance_Table.txt