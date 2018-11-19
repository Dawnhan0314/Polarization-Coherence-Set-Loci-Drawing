function ROIregion = NR_ROI_Building(pauligray)
% DESCRIPTION:����PauliRGB��ȡ��������
% OUTPUT:ROIregion:�õ��ĵ�������
% INPUT:pauligray:NL�˲�֮���PauliRGB.bmp����ת��Ϊ�Ҷ�Ӱ��

% ���Ҫ����Ľ���������
figure;imshow(pauligray);
title('���Ҫ����Ľ���������������ǵ�');
[xx,yy] = ginput(2);
% pauligraytemp = imcrop(pauligray,[xx(1),yy(1),abs(xx(1)-xx(2)),abs(yy(1)-yy(2))]);
axis([xx(1) xx(2) yy(1) yy(2)]);

% �������������Ҽ�����
title('������������');
hold on
[x,y,c]=ginput(1);
m(1)=x;
n(1)=y;
plot(x,y,'r');
k=2;
while(c==1)
    [x1,y1,c1]=ginput(1);
    if c1==1
        m(k)=x1;
        n(k)=y1;
        plot(x,y,'r');
        line([m(k-1) m(k)],[n(k-1) n(k)],'linewidth',.6,'color','r');
        k=k+1;
        c=c1;
    else
        break
    end
end
line([m(k-1) m(1)],[n(k-1) n(1)],'linewidth',.6,'color','r');
ROIregion = roipoly(pauligray,m,n);
ROIregion = double(ROIregion);
