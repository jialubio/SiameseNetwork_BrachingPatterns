% Generate training set of branches iteratively

class_ID = 1;
Nmax = 8;

for i = 1 : 100,
    ID = i;
    OptimalPattern_homogeneous(ID,class_ID,Nmax)
end
