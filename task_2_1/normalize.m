function angle_out = normalize(angle_in)
if angle_in > pi
    angle_out = angle_in - 2*pi;
elseif angle_in < - pi 
        angle_out = angle_in + 2*pi;
else
    angle_out = angle_in;
end
end

