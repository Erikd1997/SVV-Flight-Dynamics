%This program handles the calculation of Cm_alpha and Cm_delta from our
%flight test data.

%First, import the data from the excel file and convert them into SI-units:
datasheet = 'Post_Flight_Datasheet_Flight_2_DD_6_3_2018_for_test.xlsx';

ET          = xlsread(datasheet, 'C59:C65');
hp          = xlsread(datasheet, 'D59:D65').*unitsratio('meter', 'feet');
IAS         = convvel(xlsread(datasheet, 'E59:E65'), 'kts', 'm/s');
alpha       = convang(xlsread(datasheet, 'F59:F65'), 'deg', 'rad');
delta_e     = convang(xlsread(datasheet, 'G59:G65'), 'deg', 'rad');
delta_e_tr  = convang(xlsread(datasheet, 'H59:H65'), 'deg', 'rad');
F_e         = xlsread(datasheet, 'I59:I65');
FF_l        = convvel((convmass(xlsread(datasheet, 'J59:J65'), 'lbm', 'kg')), 'km/h', 'km/s');
FF_r        = convvel((convmass(xlsread(datasheet, 'K59:K65'), 'lbm', 'kg')), 'km/h', 'km/s');
F_used      = convmass(xlsread(datasheet, 'L59:L65'), 'lbm', 'kg');
TAT         = convtemp(xlsread(datasheet, 'M59:M65'), 'C', 'K');


