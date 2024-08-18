function N0 = get_NutrientConfiguration(xx, yy, rr, N0, Nmax, shape)

% GROMETRY
if shape == 1, % circle
    spot = find(rr <= 25);
    
elseif shape == 2, %ring
    spot = find(rr <= 30 & rr >= 28 ); 
    
elseif shape == 3, %cross
    spot = find((abs(xx) < 7 | abs(yy) < 7) & abs(xx) < 40 & abs(yy) < 40);
    
elseif shape == 4, %line
    spot = find(abs(xx) < 4 & abs(yy) < 40);   
    
elseif shape == 5, %hollow square
    spot = find((abs(xx - 20) < 2 | ...
        abs(yy - 20) < 2 | ...
        abs(xx + 20) < 2 | ...
        abs(yy + 20) < 2) ...
        & abs(xx) < 21 & abs(yy) < 21);
elseif shape == 6, % multiple circles
     center = [-15 0;
        15 0];
    radius = 3;
    spot = get_meshelements(xx, yy, center, radius);
    spot = find(spot == 1);
end

N0(spot) = Nmax;

end