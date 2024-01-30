function h=manifplot(manif)


rad=(manif.points.x.^2+manif.points.y.^2).^(1/2);
%-- Colormap pmin all by arclength
RGB1=[131, 195, 251]/255; % light
RGB2=[5, 52, 122]/255; % dark
R= linspace(RGB1(1),RGB2(1),100);  %// Red from 212/255 to 0
G = linspace(RGB1(2),RGB2(2),100);   %// Green from 212/255 to 0
B = linspace(RGB1(3),RGB2(3),100);  %// Blue from 1 to 170/255
RGB = [R(:), G(:), B(:)];

h=figure;
hold on
daspect([1 1 1])
view([100,30])
colormap(RGB); 
color_line3(manif.points.x,manif.points.y,manif.points.z,rad,'EdgeAlpha',1,'LineWidth',1.5);
caxis([min(rad), max(rad)]);

%-- Plane {x=y, y<0}
plane.x=[-sqrt(2)/2,0];
plane.y=[-sqrt(2)/2,0];
plane.z= [-1,1];
plane.color=[230, 178, 17]/255;
surf(repmat(plane.x,2,1), repmat(plane.y,2,1), repmat(plane.z,2,1)','FaceAlpha',0.4, 'EdgeColor',plane.color,'FaceColor',plane.color,'FaceLighting','gouraud','LineWidth',1.7)

%-- Unit circle
[xunit,yunit] = circle(0,0,1,1000);
plot3(xunit,yunit,ones(size(xunit)),'k','LineWidth',1.5)
plot3(xunit,yunit,-ones(size(xunit)),'k','LineWidth',1.5)
% 

% % %-- intersection points
idx=manif.inter.idx(manif.inter.idx<numel(manif.points.x));
plot3(manif.points.x(idx),manif.points.y(idx),manif.points.z(idx),'.','color',[57 106 160]/255,'MarkerSize',11);

%-- fixed points
plot3(manif.inf_sys.fixp.pplu.x, manif.inf_sys.fixp.pplu.y, manif.inf_sys.fixp.pplu.z,'marker','o','MarkerFaceColor',[230, 178, 17]/255,'MarkerEdgeColor',[87, 67, 6]/255,'LineWidth',1.4,'MarkerSize',6.5)
plot3(manif.inf_sys.fixp.pmin.x,manif.inf_sys.fixp.pmin.y,manif.inf_sys.fixp.pmin.z,'marker','o','MarkerFaceColor',[230, 178, 17]/255,'MarkerEdgeColor',[87, 67, 6]/255,'LineWidth',1.4,'MarkerSize',6.5)

xlabel('x')
ylabel('y')
zlabel('z')
xlim([-1.01 1.01])
ylim([-1.01 1.01])
zlim([-1.01 1.01])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = color_line3(x, y, z, c, varargin)
% color_line3 plots a 3-D "line" with c-data as color
%
%       h = color_line(x, y, z, c)
%       by default: 'LineStyle','-' and 'Marker','none'
%
%          or
%       h = color_line(x, y, z, c, mark) 
%          or
%       h = color_line(x, y, z, c, 'Property','value'...) 
%             with valid 'Property','value' pairs for a surface object
%
%  in:  x      x-data
%       y      y-data
%       z      z-data
%       c      4th dimension for colouring
%       mark   for scatter plots with no connecting line
%
% out:  h   handle of the surface object
h = surf(...
  'XData',[x(:) x(:)],...
  'YData',[y(:) y(:)],...
  'ZData',[z(:) z(:)],...
  'CData',[c(:) c(:)],...
  'FaceColor','none',...
  'EdgeColor','flat',...
  'Marker','none');
  
if nargin ==5
    switch varargin{1}
        case {'+' 'o' '*' '.' 'x' 'square' 'diamond' 'v' '^' '>' '<' 'pentagram' 'p' 'hexagram' 'h'}
            set(h,'LineStyle','none','Marker',varargin{1})
        otherwise
            error(['Invalid marker: ' varargin{1}])
    end
elseif nargin > 5
    set(h,varargin{:})
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xunit,yunit] = circle(x,y,r,n)
th = linspace(0,2*pi,n);
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
end

end