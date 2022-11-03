clear;
clc;
% x_... profile 
% y_... profile 
% V_... Vertices
% p_... point
% n_ one index on Vertices representing a point
% ns_ list of index on Vertices representing a line
% F_... Vector or matrix on Vertieces representing a face
% Parameters of gear
m_n= 2;
z= 40;
alpha_n= 20*pi/180;
x_coef=  +0;
% beta= 0*pi/180;
h_fP_coef= 1.25; % (1 + c) with c=0.25
h_aP_coef= 1.0;
rho_fP_coef= 0.25;
% s_pr= 0.0;
k= 0;
b=5;

% tip and foot radius for d_ch_flat / d_ch_down
r_a= z*m_n/2 + m_n*(h_aP_coef + x_coef + k);
r_f= z*m_n/2 - m_n*(h_fP_coef - x_coef);

% calculating gear parameters
dptch = m_n*z;
rptch = dptch/2;
dbase = dptch*cos(alpha_n);
rbase = dbase/2;

%tooth thickness eq
RY = linspace(rbase,r_a,70);
alpha_Y = acos(rbase./RY);

% Involute curve function
x1=RY.*(sin((pi/(2.*z))+(tan(alpha_n)-alpha_n)-tan(alpha_Y)+alpha_Y));
y1 = RY.*(cos((pi/(2.*z))+(tan(alpha_n)-alpha_n)-tan(alpha_Y)+alpha_Y));

% Mirroring the Curve
x1f = flip(-x1);
y1f = flip(y1);

% Topland
x2t = linspace(-x1(1,end),x1(1,end),6);
y2t = linspace(flip(y1(1,end)),y1(1,end),6);
y2tt = 0*x2t + y2t;

% chamfer parametrs
delta_a= 35*pi/180;
na_chr = x2t(end)/2;
a_ch= +na_chr;
% d_ch_flat= (2*r_f-(((r_a-r_f)*26.44)/100));
d_ch_down= (2*r_f+(((r_a-r_f)*25.33)/100));
r_e= ((((r_a-r_f)*25.33)/100)/2.3);
% a_ch = +1;
d_ch_flat = 2*r_f;
% d_ch_down = 2*r_f+2.8;
% r_e = 1.2;

% convert diameter into radius
r_ch_flat= d_ch_flat/2;
r_ch_down= d_ch_down/2;

% Calculation for fillet center
dy = cos(((pi/(2.*z))+(tan(alpha_n)-alpha_n)-tan(alpha_Y)+alpha_Y));
dx = sin(((pi/(2.*z))+(tan(alpha_n)-alpha_n)-tan(alpha_Y)+alpha_Y));

k = dy/dx;

theta = atan(-(1/k));

Ax = x1(1,1);
Ay = y1(1,1);

Cx = Ax + (r_e*cos(theta));
Cy = Ay + (r_e*sin(theta));

x3f = linspace(Ax,Cx,70);
y3f = Cy - (sqrt(r_e^2 - (x3f-Cx).^2));
y3f(1,1) = Ay;
x3ff = flip(-x3f);
y3ff = flip(y3f);
bbb = flip(y1f);
ccc = flip(x1f);
newx =(-1).*(x1f);

%deddendum circle between two adjecent teeth
phiD = 180/z;
phiDi = atand(x3f(end)/y3f(end));
rdedgap = sqrt(  (x3f(end).^2)  +  (y3f(end).^2)  );
phiDgap = (90 - (phiD+(phiD-phiDi))):1:90 - phiDi;
center = [0,0];
xgap = center(1)+(rdedgap.*cosd(phiDgap));
ygap = center(2)+(rdedgap.*sind(phiDgap));

%%%% x and y co-ordinates of complete tooth profile
x_tooth1 = [x3ff ccc x2t newx x3f xgap];
x_tooth = flip(x_tooth1);
y_tooth1 = [y3ff bbb y2t y1f y3f ygap];
y_tooth = flip(y_tooth1);

% Insert the point for d_ch_down
n_ch_down= find(r_ch_down < sqrt(x_tooth.^2 + y_tooth.^2),1,'first');
[x_tooth, y_tooth]= xy_ch_down(x_tooth, y_tooth, r_ch_down, n_ch_down);

% Add point end to the chamfer at tip diameter
n_ch_a= find(x_tooth < 0,1,'first');
x_tooth= [x_tooth(1:n_ch_a-1),a_ch,x_tooth(n_ch_a:end)];
y_tooth= [y_tooth(1:n_ch_a-1),sqrt(r_a^2 - a_ch^2),y_tooth(n_ch_a:end)];
p_ch_a = [x_tooth(n_ch_a),y_tooth(n_ch_a),0];
% Convert profile into vertices two profiles on 0 and -b height
V_tooth= [x_tooth',y_tooth', zeros(size(x_tooth')); ...
          x_tooth',y_tooth', -b*ones(size(x_tooth')) ];
% Change vertices in flat chamfer area      
V_tooth(n_ch_down:n_ch_a-1,3)= -tan(delta_a)*(V_tooth(n_ch_down:n_ch_a-1,1)-a_ch);

% Make faces for periphery
ns_top= 1:size(x_tooth');
ns_bottom= size(x_tooth')+1:size(V_tooth,1);
F_periphery= [ns_top(1:end-1)',ns_bottom(1:end-1)',ns_bottom(2:end)',ns_top(2:end)'];

% Add point ch_flat to vertices
p_ch_flat= [a_ch, sqrt(r_ch_flat^2 - a_ch^2), 0];
V_tooth= [V_tooth;p_ch_flat];
n_ch_flat= size(V_tooth,1);

% make flat chamfer face
F_ch_flat= [n_ch_down:n_ch_a,n_ch_flat];
p_ch_down= V_tooth(n_ch_down,:);

% Move and rotate the vertices so that the line from p_ch_flat to p_ch_down is in line with the y-axis!
V_tooth= V_tooth - p_ch_flat;
V_tooth= Mult_homoMat(makehgtform('yrotate',-delta_a), V_tooth);
phi_z= atan2(V_tooth(n_ch_down,1),V_tooth(n_ch_down,2));
V_tooth= Mult_homoMat(makehgtform('zrotate', phi_z), V_tooth);
V_tooth(:,3)= V_tooth(:,3)- r_e;

% rotate the vector [0,0,1] z-axis with the roation on the vertices
vec001= Mult_homoMat(makehgtform('zrotate', phi_z)* ...
                     makehgtform('yrotate',-delta_a), [0,0,1]);
                 
% calulation of the intersection between cylinder by r_e and perphery
ns_corner_down=[];
for n=n_ch_down:-1:1 
    [p_out,intersec]= Cylinder2line(r_e, V_tooth(n,:), vec001);
    if intersec==true
       V_tooth(n,:)= p_out;
       ns_corner_down= [ns_corner_down,n];
    end
end

% Move and rotate the vertices back into initial state
V_tooth(:,3)= V_tooth(:,3)+ r_e;
V_tooth= Mult_homoMat(makehgtform('zrotate', -phi_z), V_tooth);
V_tooth= Mult_homoMat(makehgtform('yrotate', delta_a), V_tooth);
V_tooth= V_tooth + p_ch_flat;

% Calculate the intesection of the cylinder r_e with x/y plane
v_ch_flat2ch_down= p_ch_flat-p_ch_down;
ns_corner_top= [1:length(ns_corner_down)]+ size(V_tooth,1);
V_tooth= [V_tooth;-V_tooth(ns_corner_down,3)/v_ch_flat2ch_down(3).*ones(size(ns_corner_down))'*v_ch_flat2ch_down + V_tooth(ns_corner_down,:)];

% make the faces
F_corner= [ns_corner_down(1:end-1)',ns_corner_top(1:end-1)',ns_corner_top(2:end)',ns_corner_down(2:end)'];
F_triangle=[ns_corner_down(end),ns_corner_top(end),ns_corner_down(end)-1];
F_bottom = ns_bottom;
F_top = [1:ns_corner_down(end)-1,ns_corner_top(end:-1:1),n_ch_a:ns_top(end)];


%tranforming the matrix for rotation
xall = (V_tooth(:,1));
yall = V_tooth(:,2);
zall = V_tooth(:,3);
xallt = xall';
yallt = yall';
zallt = zall';
%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%
r_c = 25;
deg= pi/180;
r_ce= r_c-r_e;  r_ci= 20;                                  % Start and end point of cutting edge
eta_ce= 0*deg;  eta_ci= 0*deg;                             % in polar coordinates eta / r / z
z_ce= 0;        z_ci= 0;
z_tool = -4;
delta_gear = (norm(p_ch_flat - p_ch_down)/(r_c-r_e))*z_tool/z;
v_c_comp = p_ch_flat+[delta_gear*d_ch_flat/2,0,0]-p_ch_down;
v_c_comp = v_c_comp/norm(v_c_comp);
n_ch_comp = cross(p_ch_a-(p_ch_flat+[delta_gear*d_ch_flat/2,0,0]),v_c_comp);
n_ch_comp = n_ch_comp/norm(n_ch_comp);

v_tool_0=-cross(v_c_comp, n_ch_comp)*(r_c-r_e);
p_tool_0= p_ch_flat -v_tool_0;

phi = linspace(1*delta_gear,0,30);

n_axis= n_ch_comp; % the vector of the tool axis
                 
ratio_tool_gear= -4/40;                                    % the gear ratio z_tool / z_gear

n_axis_alpha= atan2(-n_axis(2),sqrt(1-n_axis(2)^2));       % Rotaion of Z-axis (tool definition axis) into
n_axis_beta=  atan2( n_axis(1),n_axis(3));                 % rotation axis of tool in space (n_axis)

% This block calc. of tool offsetcan be replaced by a numerical calculation see below
v_tmp= makehgtform('yrotate',   -n_axis_beta)* ...         % !!! this you must adapt to your tool zero setting!
       makehgtform('xrotate',   -n_axis_alpha)* ...        % Calculation of tool offset in eta and z
       makehgtform('translate', -p_tool_0)* ...            % according to setting of p_zero_tool 
       [p_ch_flat,1]';
eta_tool_offset= atan2(v_tmp(2), v_tmp(1));   
z_tool_offset= v_tmp(3);

r_factor=linspace(1,0,10);                                % Factor to move point on cutting edge
theta_gear=   linspace(0,2*pi,300);                       % Roation of gear
[R_factor,Theta_gear]= meshgrid(r_factor,theta_gear);     % preparing surface plot with meshgrid

% here plot you gear in "zero" position
% If helicoid and gear are matching

% numerical calc. of tool offset with help of fminsearch.
[tool_offset_cal,sqr_dist] = fminsearch(@(tool_offset) sqrt_distance_XYZ(tool_offset, ...
                                                                         p_ch_flat, ...
                                                                         r_factor, theta_gear, ... 
                                                                         n_axis_alpha, n_axis_beta, ...
                                                                         p_tool_0, ...
                                                                         ratio_tool_gear, ...
                                                                         r_ce, r_ci, eta_ce, ...
                                                                         eta_ci, z_ce, z_ci), ...
                                                                         [0,0]);
disp(tool_offset_cal)
disp(sqr_dist)

% Checking results of numerical calc.
sqr_d_p_ch_flat= sqrt_distance_XYZ([eta_tool_offset,z_tool_offset], ...
                                  p_ch_flat, ...
                                  r_factor, theta_gear, ... 
                                  n_axis_alpha, n_axis_beta, ...
                                  p_tool_0, ...
                                  ratio_tool_gear, ...
                                  r_ce, r_ci, eta_ce, ...
                                  eta_ci, z_ce, z_ci);

% Finding the Z-value for a given point in x,y with help of fminsearch

p_ch_grid= [p_ch_flat(1:2),0];    % using the point p_ch_falt is just a example          
[r_factor_theta_gear,sqr_dist] = fminsearch(@(r_factor_theta_gear) ...
                                            sqrt_distance_XY (r_factor_theta_gear, ... 
                                            p_ch_grid, ...
                                            n_axis_alpha, n_axis_beta, ...
                                            p_tool_0, ...
                                            ratio_tool_gear, ...
                                            eta_tool_offset, ...                 
                                            z_tool_offset, ...
                                            r_ce, r_ci, eta_ce, ...
                                            eta_ci, z_ce, z_ci), ...
                                            [1,0]);
disp(r_factor_theta_gear)


p_ch_grid_calc= Point_on_cutting_edge(r_factor_theta_gear(1), r_factor_theta_gear(2), ...
                                      n_axis_alpha, n_axis_beta, ...
                                      p_tool_0, ...
                                      ratio_tool_gear, ...
                                      eta_tool_offset, ...
                                      z_tool_offset, ...
                                      r_ce, r_ci, eta_ce, ...
                                      eta_ci, z_ce, z_ci);
disp(p_ch_grid_calc-p_ch_flat)
%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%
% Rotational Matrix for whole gear with z teeth.
f = figure('WindowState','maximized');
pause(1);
f.Position

xlim([-16 16]);
ylim([-53 -10]);
zlim([-3 3]);
xlabel("xaxis");
ylabel("yaxis");
zlabel("zaxis");
view(4,69)

for i = 1:z
    thetar=(2*pi/z)*i;
    R = [cos(thetar) -sin(thetar) 0;sin(thetar) cos(thetar) 0;0 0 1];
    V = [xallt;yallt;zallt];
    S1 = R*V;
    xr1 = S1(1,:);
    yr1 = S1(2,:);
    zr1 = S1(3,:);
%     zr2 = 5* ( (ones(1,length(zr1))) );
%     N = patch(xr1,yr1,zr2);
    S0 = S1';
    O = patch('Vertices',S0,'Faces',F_ch_flat,'FaceColor','yellow','EdgeColor','none');
    P = patch('Vertices',S0,'Faces',F_corner,'FaceColor','blue','EdgeColor','none');
    Q = patch('Vertices',S0,'Faces',F_triangle,'FaceColor','yellow','EdgeColor','none');
    R = patch('Vertices',S0,'Faces',F_periphery,'FaceColor',[0.7,0.7,0.7],'EdgeColor','none');
    S = patch('Vertices',S0,'Faces',F_top,'FaceColor',[0.7,0.7,0.7],'EdgeColor','none');
    T = patch('Vertices',S0,'Faces',F_bottom,'FaceColor',[0.7,0.7,0.7],'EdgeColor','none');
    mgear(i,:) =[O P Q R S T];
%     plot3(S0(:,1),S0(:,2),S0(:,3));
end
    lightangle(gca, -45, 30);
    lightangle(gca, -210, 30);
    lightangle(gca, -210, 30);
hold on
r=37.1;
h=-5;
ang=0;
[a,b,c]=cylinder([0,r,r,0],100);
c([1,2],:)=0;
c([3,4],:)=h;
plot = surf(a,b,c,'FaceColor',[0.7,0.7,0.7],'EdgeColor','none');
rotate(plot, [1 1 0], ang)

for n= 1:size(R_factor,1)                                 % Loop to caculate Helicoid surface
    for m= 1:size(R_factor,2)
        p_tmp = Point_on_cutting_edge(R_factor(n,m), Theta_gear(n,m), ...
                                      n_axis_alpha, n_axis_beta, ...
                                      p_tool_0, ...
                                      ratio_tool_gear, ...
                                      eta_tool_offset, ...
                                      z_tool_offset, ...
                                      r_ce, r_ci, eta_ce, ...
                                      eta_ci, z_ce, z_ci);
        P_cx(n,m)= p_tmp(1); P_cy(n,m)= p_tmp(2); P_cz(n,m)= p_tmp(3);
    end
 
    plot3(P_cx, P_cy, P_cz);
    drawnow;
end
surf(P_cx, P_cy, P_cz)
function [x_tooth, y_tooth]= xy_ch_down(x_tooth, y_tooth, r_ch_down, n_ch_down)
% Calculates with a 4 point cubic spline the point p_ch_down
% add the point p_ch_down to the profile
   [phi_tooth,r_tooth]= cart2pol(x_tooth(n_ch_down-2:n_ch_down+1),y_tooth(n_ch_down-2:n_ch_down+1));
   phi_ch_down=spline(r_tooth, phi_tooth, r_ch_down);
   x_tooth= [x_tooth(1:n_ch_down-1),r_ch_down*cos(phi_ch_down),x_tooth(n_ch_down:end)];
   y_tooth= [y_tooth(1:n_ch_down-1),r_ch_down*sin(phi_ch_down),y_tooth(n_ch_down:end)];
end

function V= Mult_homoMat(M,V)
% Matrix multiplication on [4x4] with Vextor area [nx3] 
   V= (M * [V,ones(size(V(:,1)))]')';
   V= V(:,1:3);
end

function [p_out,intersec]= Cylinder2line(r_e, p_in, vec)
   % Calculate the insection of a line given by a point of the line
   % and the direction of the line with 
   % a cylinder with r_e around the Y-axis
   a_x= p_in(1); a_z= p_in(3);
   n_x= vec(1);  n_z= vec(3);
   tmp= (- a_x^2*n_z^2 + 2*a_x*a_z*n_x*n_z - a_z^2*n_x^2 + n_x^2*r_e^2 + n_z^2*r_e^2);
   if tmp < 0
      intersec= false; p_out=[NaN,NaN,NaN];
   else
      intersec= true;
      l= -(a_x*n_x + a_z*n_z + sqrt(tmp))/(n_x^2 + n_z^2);
      p_out= l*vec + p_in;
   end
end
 function p_c = Point_on_cutting_edge(r_factor, ...                        % [1,0] 1 outer radius, 0 inner radius
                                            theta_gear, ...                      % rotation of gear [rad]
                                            n_axis_alpha, n_axis_beta, ...       % Inclination of tool axis again z-axis first Ry(beta) then Rx(alpha)
                                            p_tool_0, ...                        % Zero point of tool
                                            ratio_tool_gear, ...                 % z_tool / z_gear
                                            eta_tool_offset, ...                 % offset of phi
                                            z_tool_offset, ...                   % offset of z
                                            r_ce, r_ci, ...                      % out, inner cutting edge point
                                            eta_ce, eta_ci, ...                  % polar coordinates  
                                            z_ce, z_ci)                          % r / eta / z
   eta_ce= eta_ce + eta_tool_offset;                                             % apply offset in eta
   eta_ci= eta_ci + eta_tool_offset;
   p_c= [r_factor*(r_ce*cos(eta_ce) - r_ci*cos(eta_ci)) + r_ci*cos(eta_ci), ...  % calculate point on cutting edge 
         r_factor*(r_ce*sin(eta_ce) - r_ci*sin(eta_ci)) + r_ci*sin(eta_ci), ...  % from polar into cart. including r_factor
         r_factor.*(z_ce - z_ci) + z_ci + z_tool_offset,1];
   p_c= makehgtform('zrotate',  -theta_gear)* ...                                % gear rotation transfered to tool
        makehgtform('translate', p_tool_0)* ...                                  % zero point of tool
        makehgtform('xrotate',   n_axis_alpha)* ...                              % apply tool axis
        makehgtform('yrotate',   n_axis_beta)* ...                              
        makehgtform('zrotate',   theta_gear/ratio_tool_gear)*p_c';               % tool roation
   p_c= p_c(1:3)';
end

function sqr_d= sqrt_distance_XYZ(tool_offset, ...
                                  p_ch_flat, ...
                                  r_factor, theta_gear, ... 
                                  n_axis_alpha, n_axis_beta, ...
                                  p_tool_0, ...
                                  ratio_tool_gear, ...
                                  r_ce, r_ci, eta_ce, ...
                                  eta_ci, z_ce, z_ci)
   
   p_c = Point_on_cutting_edge(1, 0, ...
                               n_axis_alpha, n_axis_beta, ...
                               p_tool_0, ...
                               ratio_tool_gear, ...
                               tool_offset(1), ...
                               tool_offset(2), ...
                               r_ce, r_ci, eta_ce, ...
                               eta_ci, z_ce, z_ci);
   
   sqr_d= sum((p_c - p_ch_flat).^2,'All');
end

function sqr_d= sqrt_distance_XY (r_factor_theta_gear, ... 
                                  p_ch_grid, ...
                                  n_axis_alpha, n_axis_beta, ...
                                  p_tool_0, ...
                                  ratio_tool_gear, ...
                                  eta_tool_offset, ...                 
                                  z_tool_offset, ...
                                  r_ce, r_ci, eta_ce, ...
                                  eta_ci, z_ce, z_ci)
   
   p_c = Point_on_cutting_edge(r_factor_theta_gear(1), r_factor_theta_gear(2), ...
                               n_axis_alpha, n_axis_beta, ...
                               p_tool_0, ...
                               ratio_tool_gear, ...
                               eta_tool_offset, ...                 
                               z_tool_offset, ...
                               r_ce, r_ci, eta_ce, ...
                               eta_ci, z_ce, z_ci);
   
   sqr_d= sum((p_c(1:2) - p_ch_grid(1:2)).^2,'All');
end
