function R = my_rot(theta, axis)
c = cos(theta);
s = sin(theta);

switch axis 
    
    case 'x'
        R = [1 0 0; 0 c -s; 0 s c];
    
    case 'y'
        R = [c 0 s; 0 1 0; -s 0 c];
     
    case 'z'
        R = [c -s 0; s c 0; 0 0 1];
    otherwise 
        disp("Asse non Valido");
        R = zeros(3, 3);
end