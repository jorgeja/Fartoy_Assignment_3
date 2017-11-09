function angle_out = normalize(angle_in)
for i = 1:length(angle_in)
    angle_in(i) = mod(angle_in(i),2*pi);
    if angle_in(i) > pi
        angle_out(i) = angle_in(i) - 2*pi;
    elseif angle_in(i) < - pi 
        angle_out(i) = angle_in(i) + 2*pi;
    else
        angle_out(i) = angle_in(i);
    end
end
end

