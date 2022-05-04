filename = 'data/DailyDelhiClimateTrain.csv';
M = readmatrix(filename);

% cols in weather are meantemp, humidity, wind_speed, mean_pressure
weather = M(:,2:end);
TEMP = 1;
HUMIDITY= 2;
WINDSPEED = 3;
PRESSURE = 4;
weather = rmoutliers(weather);
length(weather)

% plotting each data point
subplot(2,2,1)
plot(weather(:,TEMP))
title('Temperature')

subplot(2,2,2)
plot(weather(:,HUMIDITY))
title('Humidity')

subplot(2,2,3)
plot(weather(:,WINDSPEED))
title('Windspeed')

subplot(2,2,4)
plot(weather(:,PRESSURE))
title('Pressure')