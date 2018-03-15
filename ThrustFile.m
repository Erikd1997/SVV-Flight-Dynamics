function T = ThrustFile(c,hp_col, M_col, Mf1_col, Mf2_col, T_temp_col)
%This function obtains the thrust per engine for each moment in time
%% Firstly, the required parameters for the thrust.exe file are calculated
T_temp_ISA = c.Temp0 + c.lambda*hp_col;

Delta_T_temp_col = T_temp_col - T_temp_ISA;

% fid = fopen('matlab.dat','w');
% textarray = string(num2str(hp_col,4))+' '+string(num2str(M_col,4))+' '+string(num2str(Mf1_col,4))+' '+string(num2str(Mf2_col,4))+' '+string(num2str(Delta_T_temp_col,4));
% text = textarray(1);
% for i = 2:length(textarray)
%     text = text + '\r\n' + textarray(i);
% end
% %The matlab.dat file is written and saved
% fprintf(fid,text);
% fclose(fid);

%% Run the thrust calculations
!Thrust.exe

%% Save the thrust values
fid = fopen('Thrust.dat','r');
Result = fscanf(fid,'%f',[2 length(textarray)]);
T = Result(1,:)';
end