SEQUENCE = './sequence/';
START = 1;
% charge le nom des images de la sequence
filenames = dir([SEQUENCE '*.png']);
filenames = sort({filenames.name});
T = length(filenames);
% charge la premiere image dans 'im'
tt = START;
im = imread([SEQUENCE filenames{tt}]);
% affiche 'im'
figure;
set(gcf,'DoubleBuffer','on');
imagesc(im);