

%% Lecture de des imagesde la sequence
SEQUENCE = './seq1/';

% charge le nom des images de la sequence
filenames = dir([SEQUENCE '*.png']);
filenames = sort({filenames.name});
T = length(filenames);

img1 = imread([SEQUENCE filenames{1}]);
[Xdim, Ydim, Cdim] = size(img1);

%% Initialisation des paramètres
N = 200;
Nb = 10;
lambda = 100;
c1 = 1000;
c2 = 1000;
c3 = 10;

%% Affichage de la première image
figure
set(gcf,'DoubleBuffer','on');
imagesc(img1);

%% Selection et affichage
disp('Cliquer 4 points dans l''image pour definir la zone a suivre.');

% on recupere la zone a tracker
zone = zeros(2,4);
compteur=1;
while(compteur ~= 5)
    [x,y,button] = ginput(1);
    zone(1,compteur) = x;
    zone(2,compteur) = y;
    text(x,y,'X','Color','r');
    compteur = compteur+1;
end
newzone = zeros(2,4);
newzone(1,:) = sort(zone(1,:));
newzone(2,:) = sort(zone(2,:));

% definition de la zone a tracker
% x haut gauche, y haut gauche, largeur, hauteur
zoneAT = zeros(1,4);
zoneAT(1) = newzone(1,1);
zoneAT(2) = newzone(2,1);
zoneAT(3) = newzone(1,4)-newzone(1,1);
zoneAT(4) = newzone(2,4)-newzone(2,1);

% affichage du rectangle
rectangle('Position',zoneAT,'EdgeColor','r','LineWidth',3);

%% Calcul de l'histogramme de couleur de référence
window = imcrop(img1, zoneAT(1:4));
[tmp, Cmap] = rgb2ind(window, Nb, 'nodither');
window = rgb2ind(window, Cmap, 'nodither');
histoRef = imhist(window, Cmap);
histoRef = histoRef/sum(histoRef);

% Affichage de l'histograme de ref
%figure
%imhist(tmp, Cmap);

%% Initialisation des particules
x_0 = zoneAT(1) + zoneAT(3)/2.0
y_0 = zoneAT(2) + zoneAT(4)/2.0

sigma_0 = 0.00001;

particles = sigma_0 * randn(3,N) + repmat([x_0 ; y_0 ; 0],1,N);
particles(3,:) = ones(1,N) * 100.0;

p_w = ones(1,N);

C = [sqrt(c1) 0 0 ; 0 sqrt(c2) 0 ; 0 0 sqrt(c3)];

for t = 1:T
    img = imread([SEQUENCE filenames{t}]);
    % (a) Diffusion des particules
    for i = 1:N
        in_img = false;
        while not(in_img)
            new_particle = particles(:,i) + C*randn(3,1);
            
            x = new_particle(1);
            y = new_particle(2);
            s = new_particle(3);
            
            in_img = true;
            if floor(x - s/100.0 * zoneAT(3) / 2.0) < 1
                in_img = false;
            end
            if floor(y - s/100.0 * zoneAT(4) / 2.0) < 1
                in_img = false;
            end
            if floor(x + s/100.0 * zoneAT(3) / 2.0) > Ydim
                in_img = false;
            end  
            if floor(y + s/100.0 * zoneAT(4) / 2.0) > Xdim
                in_img = false;
            end
        end
        particles(:,i) = new_particle;
    end
    % (b) Calcul des histogrammes associées et des poids d'importance
    clf;
    imagesc(img);
    for i = 1:N
        x = particles(1,i);
        y = particles(2,i);
        s = particles(3,i);
        
        zone_p = zeros(1,4);
        zone_p(1) = floor(x - s/100.0 * zoneAT(3) / 2.0);
        zone_p(2) = floor(y - s/100.0 * zoneAT(4) / 2.0);
        zone_p(3) = floor(s/100.0 * zoneAT(3) / 2.0);
        zone_p(4) = floor(s/100.0 * zoneAT(4) / 2.0);
        
        % (d) Affichage des particules
        rectangle('Position',zone_p,'EdgeColor','b','LineWidth',1); 
        
        impart = imcrop(img, zone_p(1:4));
                 
        impart = rgb2ind(impart, Cmap, 'nodither');
        histo = imhist(impart, Cmap);
        histo = histo/sum(histo);
        
        D = sqrt( 1 - sum(sqrt(histo .* histoRef)));
        
        g = exp(-lambda * D^2);
        
        p_w(1,i) = g;                 
    end
    
    % (c) Normalisation des poids d'importance
    p_w;
    sum(p_w);
    p_w = p_w / sum(p_w);
    
    % (d) Affichage de l'estimation  
    estimation = sum(particles .* repmat(p_w, 3, 1), 2);
    x = estimation(1);
    y = estimation(2);
    s = estimation(3);

    zone_e = zeros(1,4);
    zone_e(1) = floor(x - s/100.0 * zoneAT(3) / 2.0);
    zone_e(2) = floor(y - s/100.0 * zoneAT(4) / 2.0);
    zone_e(3) = s/100.0 * zoneAT(3);
    zone_e(4) = s/100.0 * zoneAT(4);
    rectangle('Position',zone_e,'EdgeColor','r','LineWidth',2)
    drawnow;
    
    % (e) Ré-échantillonage
    new_particles = zeros(3,N);

    cumulative_p = cumsum(p_w);
    for k = 1:N
        r = rand;
        j = 1;
        while cumulative_p(j) <= r
            j = j + 1;
        end
        new_particles(:,k) = particles(:, j);  
    end

    particles = new_particles;
end







