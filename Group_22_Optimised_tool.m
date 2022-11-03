clear all;

r_collar = 7;
x_centre_s_c =  0;
y_centre_s_c = -26;
th_s_c = 0:pi/50:2*pi;
x_all_collar = r_collar * cos(th_s_c) + x_centre_s_c;
y_all_collar = r_collar * sin(th_s_c) + y_centre_s_c;
z_all_collar = ones(size(th_s_c));
l_collar = 100;
b = -(l_collar/2);
w = l_collar/2;
zz = b.*z_all_collar;
zz1 = w.*z_all_collar;

r_shaft = 35/2;
x_centre_s_c =  0;
y_centre_s_c = -26;
th_s_c = 0:pi/50:2*pi;
x_all_shaft = r_shaft * cos(th_s_c) + x_centre_s_c;
y_all_shaft = r_shaft * sin(th_s_c) + y_centre_s_c;
z_all_shaft = ones(size(th_s_c));
l_shaft = 60;
b1 = -(l_shaft/2);
w1 = l_shaft/2;
zz11 = b1.*z_all_shaft;
zz111 = w1.*z_all_shaft;
% figure(1)
% plot3(x_all_collar, y_all_collar, zz)
% grid on
% axis([-1.5  1.5    -1.5  1.5    0  10])



%%%%%%%%%%%%% INPUT FROM THE USER %%%%%%%%%%%%%%%%
d_c=60;
r_c=d_c/2;
d_min=35;
r_min=d_min/2;
d_kink=40;
r_kink=d_kink/2;
r_e=4;
phi_left=20;
phi_right = 30;
phi_kink=9;

Y=-(r_c-r_e-r_min);
Yk=-(r_c-r_e-r_kink);

%%%%%%%%%%%%% RIGHT LINE %%%%%%%%%%%%%%%%
%y1=m1*x1+c1 where:-
m1=tand(90+phi_right);
% Assuming the point (a, b) as a tangent point, hence;-
m1c=-1/m1;
%m1c- slope of radius of the circle connecting center and (a,b)
B=sqrt(r_e^2/((m1^2)+1));
A=B/m1c;
c1=-(m1*A - B);
% to find the first starting point; (C,r_min);
C=(Y-c1)/m1;
x1=linspace(C,A);
y1=(m1*x1)+c1;
% plot(x1,y1);

%%%%%%%%%%%%% LEFT LINE %%%%%%%%%%%%%%%%
%y2=m2x2+c2
m2=tand(90-phi_left);
H=sqrt(r_e^2/((m2^2)+1));
G=-sqrt(r_e^2 - H^2);
c2=H-(m2*G);
E=(Yk-c2)/m2;
x2=linspace(E,G);
y2=m2*x2+c2;
% plot(x2,y2)

%%%%%%%%%%%%% LEFT LINE  2 %%%%%%%%%%%%%%%%
%y3=m*x3+c3;
m3=tand(90-phi_kink);
c3=Yk-(m3*E);
F=(Y-c3)/m3;
x3=linspace(F,E);
y3=m3*x3+c3;
% plot(x3,y3);

%%%%%%%%%%%%% Curve %%%%%%%%%%%%%%%%
% theta=0:0.01:2*pi;
% center=[0,0];
% xc=center(1)+r_e*cos(theta);
% yc=center(2)+r_e*sin(theta);
xc=linspace(G,A);
yc=sqrt((r_e.^2)-xc.^2);
% plot(xc,yc);

%%%%%%%%%%%%% St LINE %%%%%%%%%%%%%%%%
x4=linspace(F,C);
m4=tand(0);
y4=m4*x4+Y;
% plot(x4,y4)

x1f = flip(x1);
x4f = flip(x4);
y1f = flip(y1);
%%%%%%%%%%%%% PATCH %%%%%%%%%%%%%%%%
yall=[y4 y3 y2 yc y1f];
xall=[x4f x3 x2 xc x1f];
zall= ones(1,length(xall));

% plot3(xall, yall, zall)
y3s = linspace(y4(end),y2(1));
y2s = linspace(y3s(end),((yc(1))/1.5));
ycs = yc./1.5;
y1fs = linspace(ycs(end),y4(1));
ys=[y4 y3s y2s ycs y1fs];
xs = xall./1.5;
zs = 3.*zall;
% plot3(xs,yall,zs)

figure();
xlim([-60 60]);
ylim([-70 30]);
zlim([-40 40]);
xlabel("xaxis");
ylabel("yaxis");
zlabel("zaxis");
view(161,-74)
lightangle(gca, -45, 30);
lightangle(gca, -210, 30);
lightangle(gca, -210, -60);
title('Optimised Tool Geometry....phi_right is Customer Input ### Phi_Right range : 0 - 30')
for t = 1:((length(x_all_collar))-1)
    xr = [x_all_collar(t), x_all_collar(t), x_all_collar(t+1), x_all_collar(t+1)];
    yr = [y_all_collar(t), y_all_collar(t), y_all_collar(t+1), y_all_collar(t+1)];
    zr = [zz(t), zz1(t), zz1(t+1), zz(t+1)];
    s = patch(zr,yr,xr,'cyan','EdgeColor','none');
    d(t,:) = (s);
    hold on;
end

% figure(1)
% plot3(x_all_collar, y_all_collar, zz)
% grid on
% axis([-1.5  1.5    -1.5  1.5    0  10])


for t = 1:((length(x_all_shaft))-1)
    xr = [x_all_shaft(t), x_all_shaft(t), x_all_shaft(t+1), x_all_shaft(t+1)];
    yr = [y_all_shaft(t), y_all_shaft(t), y_all_shaft(t+1), y_all_shaft(t+1)];
    zr = [zz11(t), zz111(t), zz111(t+1), zz11(t+1)];
    s = patch(zr,yr,xr,'blue','EdgeColor','none');
    e(t,:)= (s);
    hold on;
end
c1 = patch(zz, y_all_collar, x_all_collar,'cyan','EdgeColor','none');
hold on;
c2 = patch(zz1, y_all_collar, x_all_collar,'cyan','EdgeColor','none');
hold on;
c3 = patch(zz11, y_all_shaft, x_all_shaft,'blue','EdgeColor','none');
hold on;
c4 = patch(zz111, y_all_shaft, x_all_shaft,'blue','EdgeColor','none');
hold on;

for t = 1:(length(xall)-1)
    xr1 = [xall(t), xs(t), xs(t+1), xall(t+1)];
    yr1 = [yall(t), ys(t), ys(t+1), yall(t+1)];
    zr1 = [zall(t), zs(t), zs(t+1), zall(t+1)];
    s = patch(xr1,yr1,zr1,'yellow','EdgeColor','none');
    direction = [0 1 0];
    rotate(s,direction,25)
    f(t,:)=(s);
    hold on;
end
s0 = patch(xall,yall,zall,'red','EdgeColor','none');
hold on;
s1 = patch(xs,ys,zs,'red','EdgeColor','none');
hold off;
% patch(XX,YY,ZZ,'yellow');
% patch(XX,YY,zz1,'red');
direction = [0 1 0];
    rotate(s1,direction,25)
    rotate(s0,direction,25)
ii=0;
while ii < 360
      drawnow;
%        direction = [0 0 1];
%        rotate(d,direction,2)
%        rotate(s0,direction,2)
%        rotate(s1,direction,2)
%     rotate(m,direction,2)
%     origin = [-1.8,-50.99,+0.15];
%     direction = [0.5636 0.0 -0.8191];
    origin = [0, -26,  0];
    direction = [1,0,0];
    %rotate(d,direction,2,origin)
    rotate(s0,direction,-20,origin)
    rotate(s1,direction,-20,origin)
    rotate(d,direction,-20,origin)
    rotate(e,direction,-20,origin)
    rotate(f,direction,-20,origin)
    rotate(c1,direction,-20,origin)
    rotate(c2,direction,-20,origin)
    rotate(c3,direction,-20,origin)
    rotate(c4,direction,-20,origin)
 
    drawnow;
    ii = ii + 1;
end

% for i = 1:0.001:360
%      direction = [0 0 1];
%     rotate(s1,i,25)
%     rotate(s0,i,25)
%     drawnow;
%     pause(0.1)
% end
