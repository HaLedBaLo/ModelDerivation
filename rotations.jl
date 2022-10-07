function rot_x(angle)
    return [
        1 0           0;
        0 cos(angle) -sin(angle);
        0 sin(angle)  cos(angle)
    ]
end

function rot_y(angle)
    return [
         cos(angle) 0 sin(angle);
         0          1          0;
        -sin(angle) 0 cos(angle)
    ]
end

function rot_z(angle)
    return [
        cos(angle) -sin(angle) 0;
        sin(angle)  cos(angle) 0;
        0           0          1
    ]
end

function rot_3d(angle_x, angle_y, angle_z)
    return rot_z(angle_z) * rot_y(angle_y) * rot_x(angle_x)
end