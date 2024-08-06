% Open-loop transfer function parameters
A0 = 100;  % DC gain
f2 = 1e10; % Second pole frequency (Hz)
w2 = 2*pi*f2;  % Second pole frequency (rad/s)

% Phase margin values (in degrees)
PM = [60, 45, 30];

% Corresponding ωd values (rad/s)
wd_values = [3.63e8, 6.28e8, 1.088e9];

% Time vector for step response
t = 0:1e-11:5e-9;

% Plot step responses for each ωd value
figure;
hold on;

for i = 1:length(PM)
    wd = wd_values(i);
    
    % Open-loop transfer function
    s = tf('s');
    A = A0 / (1 + s/wd) / (1 + s/w2);
    
    % Closed-loop transfer function
    T = feedback(A, 1);
    
    % Plot step response
    step(T, t);
end

legend('PM = 60°', 'PM = 45°', 'PM = 30°');
xlabel('Time (s)');
ylabel('Amplitude');
title('Step Response for Different Phase Margins');

% Set axis limits to zoom in
xlim([0, 5e-9]);
ylim([0, 2]);

grid on;
hold off;