function [ROIrow,ROIcol] = NR_ROI_Point_Building(pauligray)
% DESCRIPTION:利用PauliRGB提取叠掩区域
% OUTPUT:ROIregion:得到的叠掩区域
% INPUT:pauligray:NL滤波之后的PauliRGB.bmp，并转化为灰度影像

% 框出要计算的建筑物区域
figure;imagesc(pauligray);
% title('框出要计算的建筑物区域的两个角点');
% [xx,yy] = ginput(2);
% % pauligraytemp = imcrop(pauligray,[xx(1),yy(1),abs(xx(1)-xx(2)),abs(yy(1)-yy(2))]);
% axis([xx(1) xx(2) yy(1) yy(2)]);

% 勾画感兴趣点，右键结束
title('勾画兴趣点');
hold on
[ROIrow,ROIcol,c]=ginput(1);
ROIrow = round(ROIrow);
ROIcol = round(ROIcol);
% k=2;
% while(c==1)
%     [x1,y1,c1]=ginput(1);
%     if c1==1
%         m(k)=x1;
%         n(k)=y1;
%         plot(x,y,'r');
%         line([m(k-1) m(k)],[n(k-1) n(k)],'linewidth',.6,'color','r');
%         k=k+1;
%         c=c1;
%     else
%         break
%     end
% end
% line([m(k-1) m(1)],[n(k-1) n(1)],'linewidth',.6,'color','r');
