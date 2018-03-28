function [F_used1, F_used2, F_used3, payload, Fuel0] = FuelReader(datasheet)
F_used1 = xlsread(datasheet, 'I28:I33');
F_used2 = xlsread(datasheet, 'L59:L65');
F_used3 = xlsread(datasheet, 'L75:L76');

payload = xlsread(datasheet, 'H8:H16');
Fuel0 = xlsread(datasheet, 'D18:D18');