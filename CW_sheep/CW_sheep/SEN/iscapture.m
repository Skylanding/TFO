function captured = isCaptured(particle_position, direction, sensor_position, organ_radius)
    % 计算粒子放射方向的单位向量
    dir_vector = [cos(direction), sin(direction)];
    
    % 计算传感器位置的单位向量
    sensor_vector = sensor_position / norm(sensor_position);
    
    % 计算粒子放射方向与传感器位置向量之间的角度
    angle = acos(dot(dir_vector, sensor_vector));
    
    % 检查这个角度是否小于传感器覆盖的角度范围
    % 假设传感器覆盖的角度范围是均匀分布的
    sensor_angle_span = 2 * pi / num_sensors;
    
    captured = angle < sensor_angle_span / 2;
end
