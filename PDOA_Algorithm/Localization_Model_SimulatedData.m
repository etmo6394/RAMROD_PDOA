clear all; close all; clc;

%% Parameters
n_samples = 200;%number of samples along path

wavelength = 1.7241e-10;%Wavelength [m]
P_t = 0.6;%Power [W]
alpha = 4;%Path loss coeff
G_t = 0.8;%Transmitter Gain
G_r = 0.8;%Receiver Gain
L = 1;%Additional losses
path_error = 100;%Error in x and y positions [m]

% Random RFI Source Location
x_source = randi(3000,1);
y_source = randi(3000,1);
middle = 1500;

%Radius of influence RFI
r = 500;
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);

% Define UAS Path
x = linspace(0,3000,n_samples);
y = 1500+1500*sin(0.02*x) + (-path_error/2 + (path_error).*rand(1,n_samples));
z = linspace(40,40,n_samples);
dist_path = sqrt((x_source-x).^2 + (y_source-y).^2);%distance b/t source and path [m]
P_actual_pathW = P_t*G_t*G_r*(wavelength^2)./((4*pi)^2.*dist_path.^alpha*L);%Power along path [m]
P_actual_path = 10*log10(1000*P_actual_pathW);%Power along path [dBm]

% Take out data not in 500m radius
j = 1;
for i = 1:length(x)
	if sqrt((x(i)-x_source)^2 + (y(i)-y_source)^2) <= 500
        %Path data within radius of influence
        x_data(j) = x(i);
        y_data(j) = y(i);
        
        %Path data with random error
        x_spoof(j) = x_data(j) + (-path_error/2 + (path_error).*rand(1));
        y_spoof(j) = y_data(j) + (-path_error/2 + (path_error).*rand(1));
        
        P_actual_path_data(j) = P_actual_path(i);%power along path [dBm]
        j = j + 1;
    end
end
z_data = linspace(40,40,length(x_data));


% Simulated ARBITRARY Grid to Create Objective Function
x_coords = x;
y_coords = x;
[X,Y] = meshgrid(x_coords,y_coords);

%Initialize Matrices
dP_m = zeros(length(x));
dP_a = zeros(length(x));
Q = zeros(length(x_coords));
P_actual = zeros(length(x_coords));

tic
for i = 1:length(x_coords)
    for j = 1:length(y_coords)
        for l = 1:length(x_data)
            for k = 1:l
                %Calculate CALCULATED power difference b/t every set of points
                %uses location w/ error b/c the calculated data occurs
                %where we THINK we are
            	dist_l = sqrt((x_coords(i)-x_spoof(l)).^2 + (y_coords(j)-y_spoof(l)).^2);
                dist_k = sqrt((x_coords(i)-x_spoof(k)).^2 + (y_coords(j)-y_spoof(k)).^2);
                dP_m(l,k) = 5*alpha*log10(dist_l^2/dist_k^2);

                
                % ^^^^^
                % Separate Things
                % (there is no down arrow but imagine there's some here)
                
                
                %Calculate ACTUAL power difference b/t every set of points
                %****This will be replaced with actual data*****
                %The source location is used in conjunction with the path
                %loss equation to simulate power data (what our instruments will read)
                %k and l are measurement locations
                %does not use spoofed data b/c the actual power
                %measurements occur where we ACTUALLY are, not where we
                %THINK we are
                dk = sqrt((x_source-x_data(k)).^2 + (y_source-y_data(k)).^2);
                P_actual_Wk = P_t*G_t*G_r*(wavelength^2)./((4*pi)^2.*dk.^alpha*L);
                P_actual_k = 10*log10(1000*P_actual_Wk);
                dl = sqrt((x_source-x_data(l)).^2 + (y_source-y_data(l)).^2);
                P_actual_Wl = P_t*G_t*G_r*(wavelength^2)./((4*pi)^2.*dl.^alpha*L);
                P_actual_l = 10*log10(1000*P_actual_Wl);
                dP_a(l,k) = P_actual_k - P_actual_l;
                
                %Calculate objctive function - probability for each
                %arbitrary grid point that RFI is located there
                Q(i,j) = Q(i,j) + (dP_a(l,k)-dP_m(l,k)).^2;
            end
        end
    end
    disp(i)
end

%Find point minimizing Q
[minQ,I] = min(Q(:));
[XI, YI] = find(Q == minQ);
toc

limits = [0,3000, 0,3000];
%distance b/t estimate and source
distance = sqrt((x_coords(XI) - x_source)^2 + (y_coords(YI) - y_source)^2);
%% Plotting

figure(1)
surf(Y,X,Q)
axis(limits)
title('Objective Function Q(x,y)')
xlabel('x [m]')
ylabel('y [m]')

figure(2)
grad = P_actual_path_data;
scatter(x_source, y_source, 'ro', 'LineWidth', 2)
hold on
scatter(x_coords(XI), y_coords(YI), 'c*', 'LineWidth', 2);
h = cline(x_spoof,y_spoof,z_data,grad,'parula');
plot(x_source+xp,y_source+yp,'r:')
title(sprintf('UAS Path (%i samples) (Error: %.2f m)',n_samples,distance))
legend('Emitter Location','Estimate')
xlabel('X Distance (km)');
ylabel('Y Distance (km)');
axis(limits)



fprintf('Distance: %f m\n',distance)
