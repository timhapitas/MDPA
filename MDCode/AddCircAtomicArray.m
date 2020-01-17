function AddCircAtomicArray(x_rad, y_rad, X0, Y0, VX0, VY0, InitDist, Temp, Type)
global C
global x y AtomSpacing
global nAtoms
global AtomType Vx Vy Mass0 Mass1

if Type == 0
    Mass = Mass0;
else
    Mass = Mass1;
end

L = (2*x_rad - 1) * AtomSpacing;
W = (2*y_rad - 1) * AtomSpacing;

xp(1, :) = linspace(-L/2, L/2, 2*x_rad);
yp(1, :) = linspace(-W/2, W/2, 2*y_rad);

numAtoms = 0;
for i = 1:2*x_rad
    for j = 1:2*y_rad
        if ((X0 - xp(i))/x_rad)^2 + ((Y0 - yp(j))/y_rad)^2 <= (AtomSpacing)^2
            numAtoms = numAtoms+1;
            x(nAtoms + numAtoms) = xp(i);
            y(nAtoms  + numAtoms) = yp(j);
        else
            i
            j
        end
    end
end


x(nAtoms + 1:nAtoms + numAtoms) = x(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + X0;
y(nAtoms + 1:nAtoms + numAtoms) = y(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + Y0;

AtomType(nAtoms + 1:nAtoms + numAtoms) = Type;

if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;
    Vy(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb * Temp / Mass);

    Vx(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
    Vy(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
end

Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms)) + VX0;
Vy(nAtoms + 1:nAtoms + numAtoms) = Vy(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vy(nAtoms + 1:nAtoms + numAtoms)) + VY0;

nAtoms = nAtoms + numAtoms;

end
