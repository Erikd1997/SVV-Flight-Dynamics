function cg_fuel = Fuel_cg_as_function_of_fuel_mass(fuel)
fuel_list = 100:100:2800;
moment_list = [298.16,591.18,879.08,1165.42,1448.40,1732.53,2014.80,2298.84,2581.92,2866.30,3150.18,3434.52,3718.52,4003.23,4287.76,4572.24,4856.56,5141.16,5425.64,5709.90,5994.04,6278.47,6562.82,6846.96,7131.00,7415.33,7699.60,7984.34,8269.06,8554.05];

cg_fuel = zeros(1,length(fuel));

for i = 1:length(fuel)
    Fuel_left = fuel_list(fuel_list<fuel(i));
    Fuel_lower = Fuel_left(end);
    Moment_lower = moment_list(find(fuel_list == Fuel_lower));
    
    Fuel_right = fuel_list(fuel_list>fuel(i));
    Fuel_higher = Fuel_right(1);
    Moment_higher = moment_list(find(fuel_list == Fuel_higher));
    
    Moment_fuel = (fuel(i)-Fuel_lower)/(Fuel_higher-Fuel_lower)*(Moment_higher-Moment_lower) + ...
        Moment_lower;
    
    cg_fuel(i) = Moment_fuel*100/fuel(i);
end
end