function [XX,YY]=splitThat(XY)

XX=XY(1:end/2,:);
YY=XY(end/2+1:end,:);

return;