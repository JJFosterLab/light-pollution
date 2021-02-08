clear
addpath('E:\sprograms\elf');
elf_paths;

load('F:\Newly recalced\m1b Occupied sunspot - anterior view D810\filt\scene001_filt.mat')
I1 = im_filt_HDR{1};
load('F:\Newly recalced\m2b Occupied sunspot - posterior view 810\filt\scene001_filt.mat')
I2 = im_filt_HDR{1};

I = [I1 I2];
I = I(:, 1:end-1, :);
In = elf_io_correctdng(I, '', 'maxval', max(I1(:)));

figure(5); clf; hold on;
[A, E] = meshgrid(-180:.1:180, 90:-.1:-90);
[X, Y, Z] = sph2cart(deg2rad(A), deg2rad(E), 1.01*ones(size(A)));
view(60, 5);
h = surf(X, Y, Z, In, 'edgecolor', 'none'); 
axis equal
shading interp
axis off

% [X1, Y1, Z1] = sphere(1000);
% surf(X1, Y1, Z1, 'edgecolor', 'none', 'facecolor', [.7 .7 .7]);

%save('F:\Newly recalced\eyemap.mat', 'eye_az', 'eye_el')
load('F:\Newly recalced\eyemap.mat', 'eye_az', 'eye_el');
eye_az = [eye_az; eye_az(1)];
eye_el = [eye_el; eye_el(1)];

[eyex, eyey, eyez] = sph2cart(deg2rad(eye_az-90), deg2rad(eye_el), 1.02*ones(size(eye_az)));
h = plot3(eyex, eyey, eyez, 'k', 'linewidth', 4)

[eyex, eyey, eyez] = sph2cart(deg2rad(-eye_az-90), deg2rad(eye_el), 1.02*ones(size(eye_az)));
h2 = plot3(eyex, eyey, eyez, 'k', 'linewidth', 4)

view(0, 0);
rotate(h, [1 0 0], -45);
rotate(h2, [1 0 0], -45);

quiver3(0, 3, 0, 0, -6, 0, 'k', 'linewidth', 5)
view(30, 0);


pdfsave(5, 'C:\Users\Jochen\Desktop\almut1.pdf');
view(210, 0);
pdfsave(5, 'C:\Users\Jochen\Desktop\almut2.pdf');