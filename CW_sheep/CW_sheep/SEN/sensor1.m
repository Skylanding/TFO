clear;

% 定义器官区域为一个圆形，半径为1
organ_radius = 1;

% 定义传感器数量和分布
num_sensors = 8; % 假设有8个传感器
angles = linspace(0, 2*pi, num_sensors + 1);
angles(end) = []; % 去掉重复的最后一个角度
sensor_positions = organ_radius * [cos(angles); sin(angles)]';

% 模拟放射源，这里假设放射源在圆心，产生粒子数为1000
num_particles = 1000;
particle_positions = repmat([0, 0], num_particles, 1); % 所有粒子从圆心放射

% 模拟粒子的放射方向
particle_directions = rand(num_particles, 1) * 2 * pi;

% 初始化传感器捕获的粒子数
sensor_counts = zeros(num_sensors, 1);

% 模拟粒子的运动，判断是否被传感器捕获
for i = 1:num_particles
    direction = particle_directions(i);
    % 粒子的路径可以用极坐标表示 (r, theta)
    % 这里我们检查路径与传感器位置的交点
    for j = 1:num_sensors
        % 计算粒子方向与传感器位置的交点
        sensor_position = sensor_positions(j, :);
        % 检查粒子是否经过传感器位置
        if isCaptured(particle_positions(i,:), direction, sensor_position, organ_radius)
            sensor_counts(j) = sensor_counts(j) + 1;
        end
    end
end

% 输出每个传感器捕获的粒子数
disp('每个传感器捕获的粒子数：');
disp(sensor_counts);

% 定义粒子捕获函数
function captured = isCaptured(particle_position, direction, sensor_position, organ_radius)
    % 这里我们需要判断直线与圆的交点是否在传感器位置上
    % 由于所有粒子都从圆心放射，所以我们只需检查方向
    sensor_angle = atan2(sensor_position(2), sensor_position(1));
    captured = abs(wrapToPi(direction - sensor_angle)) < pi / num_sensors;
end
