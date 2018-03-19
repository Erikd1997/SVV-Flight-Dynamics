function [ET, hp, IAS, alpha, FFl, FFr, F_used, TAT, Weight] = ExcelReader_S1(datasheet, BEM)
Payload_Weight = sum(xlsread(datasheet, 'H8:H16'));
Block_Fuel     = xlsread(datasheet, 'D18:D18');
ET             = xlsread(datasheet, 'C28:C33');
hp             = xlsread(datasheet, 'D28:D33')*0.3048;
IAS            = xlsread(datasheet, 'E28:E33')*0.514444444;
alpha          = deg2rad(xlsread(datasheet, 'F28:F33'));
FFl            = xlsread(datasheet, 'G28:G33')*0.45359237/3600;
FFr            = xlsread(datasheet, 'H28:H33')*0.45359237/3600;
F_used         = xlsread(datasheet, 'I28:I33')*0.45359237;
TAT            = xlsread(datasheet, 'J28:J33') + 273.15;
Weight         = Payload_Weight + Block_Fuel + BEM - F_used;