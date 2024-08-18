ID = getenv('SLURM_ARRAY_TASK_ID'); 
ID = str2num(ID); 

class_ID = 5;
Nmax = 6;

% Given the mapping between conditions and patterns
% compute the pattern under a heterogeneous condition
dir_output = strcat('output/',num2str(class_ID));
mkdir(dir_output);

% --------------------------------------------------------
% Parameters
% config: type of initial configuration
% 1 - single spot
% 2 - multiple spot(plaintext -> binary)
%     need to choose number of character, default: 1 char
config = 1;

% shape: nutrient field shape
% 1 - circle
% 2 - ring
% 3 - cross
% 4 - line
% 5 - hollow square
% 6 - multiple circles
% 7 - linear
% 8 - uniform
shape = 8;

% --------------------------------------------------------
% Grid
L      = 90;
totalt = 100;

dt = 0.02;
nt = totalt / dt;
nx = 451; % 451;
ny = nx;
dx = L / (nx - 1); dy = dx;
x  = linspace(-L/2, L/2, nx);
y  = linspace(-L/2, L/2, ny);
[xx, yy] = meshgrid(x, y);
rr = sqrt(xx .^ 2 + yy .^ 2);

% Model parameters
gama = 7.5;

bN = 160;
DN = 9;
aC = 0.8 * 1.5;
KN = 0.8;
Cm = 0.05;

noiseamp = 0 * pi;

% Initial condition
r0 = 4;     % radius of initial cell seeding
C0 = 1.6;   % initial cell density

% Initialization
P  = zeros(nx, ny);      % Pattern
C  = zeros(nx, ny);      % Cell density
N0 = zeros(nx,ny);       % Initial nutrient

% Define nutrient field N
Nmin = 2;
% Nmax = 6;
if shape == 7, % linear
    N = (Nmax - Nmin) * xx / L + (Nmax + Nmin) / 2;
elseif shape == 8, % uniform
    N = ones(nx,ny) * Nmax;
else, % geometry
    N0 = get_NutrientConfiguration(xx, yy, rr, N0, Nmax, shape);
    nt_diffusion = 15; % change to control gradient
    N = get_NutrientField(xx,yy,Nmax, N0,nt_diffusion);
    N = N * (Nmax-Nmin);
    N = N + Nmin;
end

% Get optimal width and density given N
[Wmat, Dmat] = param2pattern_rand(N);
ntips0 = round(2 * pi * r0 * Dmat((nx+1)/2,(nx+1)/2));
ntips0 = max(ntips0, 2); % initial branch number

% --------------------------------------------------------
% INITIAL PATTERN
% center:       innoculum location
% num_seeding:  number of innoculum
if config == 1, % single spot
    center = [0,0];
    num_seeding = 1;
elseif config ==2, % multiple spots
    load([pwd '/char_list.mat']);
    load([pwd '/binary_list.mat']);
    binary_list = fliplr(binary_list);
    plaintext_length = 1; % number of char in the plaintext
    [plaintext_binary, plaintext] = RandPlainText_Generator(plaintext_length,char_list,binary_list);
    fprintf(strcat('The plaintext is: ', plaintext))
    num_seeding = sum(sum(plaintext_binary));
    plaintext_binary_base = ones(plaintext_length, 6);
    xcenter = linspace(-25,25,6);
    if plaintext_length == 1,
        ycenter = 0;
    else
        ycenter = linspace(10,-10,plaintext_length);
    end
    [a b]   = size(plaintext_binary_base);
    center  = zeros(num_seeding,2);
    idx = 1;
    for i = 1:b,
        for j = 1:a,
            if(plaintext_binary(j,i) == 1),
                center(idx,1) = xcenter(i);
                center(idx,2) = ycenter(j);
                idx = idx +1;
            end
        end
    end
end


P = get_meshelements(xx, yy, center, r0);
C(P == 1) = C0 / (sum(P(:)) * dx * dy);
C_pre = C;

% Define branch
% branch domain
BranchDomain = cell(ntips0 * num_seeding, 1);
for k = 1 : ntips0 * num_seeding,
    BranchDomain{k} = C > 0;
end

% branch angle
dE    = zeros(ntips0, 1);
rng('shuffle'); rngState = rng;
seed = rngState.Seed + uint32(feature('getpid')); rng(seed);
theta0 = pi/2 + rand;
theta = linspace(theta0, 2*pi+theta0, ntips0 + 1)';

theta = theta(1 : ntips0);
theta = repmat(theta,num_seeding,1);
delta = linspace(-1, 1, 201) * pi;

% branch center
Tipx = zeros(ntips0 * num_seeding, 1); % initial branch center x
Tipy = zeros(ntips0 * num_seeding, 1); % initial branch center y
for i = 1: ntips0 * num_seeding,
    j = ceil(i/ntips0);
    Tipx(i) = center(j,1);
    Tipy(i) = center(j,2);
end

[MatV1N,MatV2N,MatU1N,MatU2N] = Branching_diffusion(dx,dy,nx,ny,dt,DN);

% --------------------------------------------------------
% GROWTH
figure
for i = 0 : 6000,
    % -------------------------------------
    % Nutrient distribution and cell growth
    fN = N ./ (N + KN) .* Cm ./ (C + Cm) .* C;
    dN = - bN * fN;
    N  = N + dN * dt;
    NV = MatV1N \ (N * MatU1N);
    N  = (MatV2N * NV) / MatU2N;
    
    dC = aC * fN;
    C  = C + dC * dt;
    
    
    % -------------------------------------
    % Branch extension and bifurcation
    if mod(i, 0.5/dt) == 0,% plot every 25 steps
        
        Width    = interp2(xx, yy, Wmat, Tipx, Tipy); % local branch width at each tip
        dBiomass = (C - C_pre) * dx * dy;             % totoal biomass change
        
        % compute the amount of biomass accumulation in each branch
        BranchDomainSum = cat(3, BranchDomain{:});
        BranchDomainSum = sum(BranchDomainSum, 3);
        ntips = length(Tipx);
        for k = 1 : ntips
            branchfract = 1 ./ (BranchDomainSum .* BranchDomain{k});
            branchfract(isinf(branchfract)) = 0;
            dE(k) = sum(sum(dBiomass .* sparse(branchfract))); % energy change
        end
        
        % extension rate of each branch
        dl = gama * dE ./ Width;
        
        noise_max = 1.3;
        noise_min = 0.7;
        rng('shuffle'); rngState = rng;
        seed = rngState.Seed + uint32(feature('getpid')); rng(seed);
        noise_dl  = noise_min+rand*(noise_max-noise_min);
        dl = dl * noise_dl; % add noise
        if i == 0,
            dl = 2;
        end
        
        % Update bifurcation
        Density  = interp2(xx, yy, Dmat, Tipx, Tipy);
        R        = 1 ./ Density;
        TipxNew  = Tipx;
        TipyNew  = Tipy;
        thetaNew = theta;
        dlNew    = dl;
        BranchDomainNew = BranchDomain;
        for k = 1 : ntips,
            dist2othertips = sqrt((TipxNew - Tipx(k)) .^ 2 + (TipyNew - Tipy(k)) .^ 2);
            dist2othertips = sort(dist2othertips);
            if dist2othertips(2) > R(k)
                TipxNew    = [TipxNew; Tipx(k) + dl(k) * sin(theta(k) + 0.5 * pi)];
                TipyNew    = [TipyNew; Tipy(k) + dl(k) * cos(theta(k) + 0.5 * pi)];
                TipxNew(k) = TipxNew(k) + dl(k) * sin(theta(k) - 0.5 * pi);
                TipyNew(k) = TipyNew(k) + dl(k) * cos(theta(k) - 0.5 * pi);
                dlNew      = [dlNew; dl(k) / 2];
                dlNew(k)   = dl(k) / 2;
                thetaNew   = [thetaNew; theta(k)];
                BranchDomainNew{end+1} = BranchDomain{k};
            end
        end
        Tipx  = TipxNew;
        Tipy  = TipyNew;
        theta = thetaNew;
        dl    = dlNew;
        ntips = length(Tipx);
        BranchDomain = BranchDomainNew;
        
        % Determine branch extension angle
        Tipx_pre = Tipx;
        Tipy_pre = Tipy;
        if i == 0,
            Tipx = Tipx + dl .* sin(theta);
            Tipy = Tipy + dl .* cos(theta);
        else
            thetaO = ones(ntips, 1) * delta;
            TipxO  = Tipx + dl .* sin(thetaO);
            TipyO  = Tipy + dl .* cos(thetaO);
            NO     = interp2(xx, yy, N, TipxO, TipyO);
            [~, ind] = max(NO, [], 2); % find the direction with maximum nutrient
            
            for k = 1 : ntips,
                Tipx(k)  = TipxO(k, ind(k));
                Tipy(k)  = TipyO(k, ind(k));
                theta(k) = thetaO(k, ind(k));
            end
        end
        
        % Stop growing when approaching edges
        ind       = sqrt(Tipx.^2 + Tipy.^2) > 0.95 * L/2;
        Tipx(ind) = Tipx_pre(ind);
        Tipy(ind) = Tipy_pre(ind);
        
        % Fill the width of the branches
        Width = interp2(xx, yy, Wmat, Tipx, Tipy);
        for k = 1 : ntips,
            d = sqrt((Tipx(k) - xx) .^ 2 + (Tipy(k) - yy) .^ 2);
            P(d <= Width(k)/2) = 1;
            BranchDomain{k} = BranchDomain{k} | (d <= Width(k)/2);
        end
        C(P == 1) = sum(C(:)) / sum(P(:));
        C_pre = C;
        
        % -------------------------------------
        % Plot
        %         clf; ind = 1 : 2 : nx;
        %         timestep = strcat(' at t = ', num2str(i));
        %         % plot cell and center of the tips
        %         subplot 121; hold on;
        %         pcolor(xx(ind, ind), yy(ind, ind), C(ind, ind)); shading interp; axis equal;
        %         axis([-L/2 L/2 -L/2 L/2]); colormap('gray');title(strcat('Cell', timestep));
        %         set(gca,'YTick',[], 'XTick',[]);colorbar
        %         plot(Tipx, Tipy, '.', 'markersize', 5);hold off
        %         % plot nutrient
        %         subplot 122
        %         pcolor(xx(ind, ind), yy(ind, ind), N(ind, ind)); shading interp; axis equal;
        %         axis([-L/2 L/2 -L/2 L/2]); set(gca,'YTick',[], 'XTick',[]); title(strcat('Nutrient', timestep));
        %         colormap('parula'); caxis([0 Nmax]); colorbar; drawnow
    end
    
    % Save initial pattern
    %     if i == 25,
    %         filename = [pwd '/output/Initial_' num2str(ID) '.fig'];
    %         savefig(filename);
    %     end
    
end

% Save final pattern
IM = C/max(max(C));
filename2 = [pwd '/' dir_output '/Final_Class_' num2str(class_ID) '_' num2str(ID) '.jpg'];
imwrite(IM, filename2, 'jpg');


