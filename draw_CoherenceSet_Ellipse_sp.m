function [VolCoherence,GroundCoherence,fai0,ellipse]=draw_CoherenceSet_Ellipse_sp(T6,ypdh,ypdl,varargin)

%T6:PolInSAR的T6矩阵
%varargin：（1）要是两个参数，即m,n，则进行单点PolInSAR相干集在复单位圆上椭圆绘制
%          （2）要是无参数，则进行整个影像PolInSAR相干集在复单位圆上椭圆计算
%pdh:PDHigh通道复相关系数，进行判断地表相位，和ymax与ymin
%ellipse：相干集扁率


[row,col]=size(ypdh);
ellipse=zeros(row,col);

%A 整个影像相干集绘制
if(numel(varargin)==0)
    
    fai0=zeros(row,col);
    Rmax=zeros(row,col);
    Rmin=zeros(row,col);
    
    
    VolCoherence=zeros(row,col);
    GroundCoherence=zeros(row,col);


%计算相干区域的长轴，接着利用长轴与复平面单位圆交点确定地表相位
   for m=1:row
       m
       for n=1:col
        
        
        
        %@1 构建每个像素的T6矩阵
        T11=[T6(1,1,m,n),T6(1,2,m,n),T6(1,3,m,n);T6(2,1,m,n),T6(2,2,m,n),T6(2,3,m,n);T6(3,1,m,n),T6(3,2,m,n),T6(3,3,m,n)];
        
        T22=[T6(4,4,m,n),T6(4,5,m,n),T6(4,6,m,n);T6(5,4,m,n),T6(5,5,m,n),T6(5,6,m,n);T6(6,4,m,n),T6(6,5,m,n),T6(6,6,m,n)];
        
        V12=[T6(1,4,m,n),T6(1,5,m,n),T6(1,6,m,n);T6(2,4,m,n),T6(2,5,m,n),T6(2,6,m,n);T6(3,4,m,n),T6(3,5,m,n),T6(3,6,m,n)];
        
        
        dist=zeros(720,1);%存储距离
        y=zeros(720,2);
        yy=zeros(720,2);
        
        %@2 相干区域采样间隔
        for k=1:720
            th=pi/720*k;
            T=(T11+T22)/2;
            A=(exp(th*1i)*V12+exp(-th*1i)*V12')/2;
            W=-inv(T)*A;
            [d dd]=eig(W);%d是特征向量，dd是特征值
            
            %@2计算Re(y*exp(faio*i))的极大小值，找出对应的特征向量w1,w2,得到相干区域边界点
            re(1)=d(:,1)'*A*d(:,1)/(d(:,1)'*T*d(:,1))*exp(th*1i);
            re(2)=d(:,2)'*A*d(:,2)/(d(:,2)'*T*d(:,2))*exp(th*1i);
            re(3)=d(:,3)'*A*d(:,3)/(d(:,3)'*T*d(:,3))*exp(th*1i);
            re=real(re);

            %判断极大值极小值对应的特征值点，计算两点之间的距离记录下
            re_max=find(re==max(re));
            re_min=find(re==min(re));
            y1=d(:,re_max)'*V12*d(:,re_max)/(d(:,re_max)'*T*d(:,re_max));
            y2=d(:,re_min)'*V12*d(:,re_min)/(d(:,re_min)'*T*d(:,re_min));
             
            dist(k)=abs(y1-y2);
            y(k,1)=angle(y1);y(k,2)=abs(y1);
            yy(k,1)=angle(y2);yy(k,2)=abs(y2);
            
        end
        
        
         
        %@3 计算最大距离对应长轴的两交点，接着计算长轴直线的参数斜距M和截距C，y=Mx+C
        ellipse(m,n)=sqrt(1-min(dist)*min(dist)/(max(dist)*max(dist)));%ellipse=sqrt(1-b*b/(a*a))
        
        num=find(dist==max(dist));
        thth=pi/720*num;
        
        T=(T11+T22)/2;
        A=(exp(thth*1i)*V12+exp(-thth*1i)*V12')/2;
        W=-pinv(T)*A;
        [d dd]=eig(W);%d是特征向量，dd是特征值
        
        re(1)=d(:,1)'*A*d(:,1)/(d(:,1)'*T*d(:,1))*exp(th*1i);
        re(2)=d(:,2)'*A*d(:,2)/(d(:,2)'*T*d(:,2))*exp(th*1i);
        re(3)=d(:,3)'*A*d(:,3)/(d(:,3)'*T*d(:,3))*exp(th*1i);
        re=real(re);   

        re_max=find(re==max(re));
        re_min=find(re==min(re));
        Rmax(m,n)=d(:,re_max)'*V12*d(:,re_max)/(d(:,re_max)'*T*d(:,re_max));
        Rmin(m,n)=d(:,re_min)'*V12*d(:,re_min)/(d(:,re_min)'*T*d(:,re_min));
        
        if(abs(Rmax(m,n)-ypdh(m,n))<abs(Rmin(m,n)-ypdh(m,n)))
            VolCoherence(m,n)=Rmax(m,n);
            GroundCoherence(m,n)=Rmin(m,n);
        else
            VolCoherence(m,n)=Rmin(m,n);
            GroundCoherence(m,n)=Rmax(m,n);
        end
       
        
        b=tan(angle(Rmax(m,n)-Rmin(m,n)));
        a=imag(Rmax(m,n))-(imag(Rmax(m,n))-imag(Rmin(m,n)))/(real(Rmax(m,n))-real(Rmin(m,n)))*real(Rmax(m,n));
        
        
        %@4 计算相干区域长轴与复平面单位圆的交点，并判断地表相位
        %计算拟合直线y=bx+a与单位圆x^2+y^2=1的两个交点
        delta=4*a^2*b^2-4*(1+b^2)*(a^2-1);
        fai1_x=(-2*a*b+sqrt(delta))/(2*(1+b^2));
        fai1_y=b*fai1_x+a;
        fai1=fai1_x+fai1_y*1i;
        fai2_x=(-2*a*b-sqrt(delta))/(2*(1+b^2));
        fai2_y=b*fai2_x+a;
        fai2=fai2_x+fai2_y*1i;
        %将某一交点与HV（体散射占优，代表植被相位）距离大于离HH-VV（二面角散射占优，代表地表相位）的距离比较，为判断地表相位点
        dist_pdh=abs(fai1-ypdh(m,n));
        dist_pdl=abs(fai1-ypdl(m,n));
%           angle1=angle(fai1);
%           angle2=angle(fai2);
        if(dist_pdh>dist_pdl)
            fai0(m,n)=angle(fai1);
        else
            fai0(m,n)=angle(fai2);
        end
        
        
       end
    end    
    
    
    
    
%B 单点相干集绘制
else if(numel(varargin)==2)
        
        
        %获取单点坐标
        m=varargin{1};
        n=varargin{2};
        
        
        
        %@1 构建每个像素的T6矩阵
        T11=[T6(1,1,m,n),T6(1,2,m,n),T6(1,3,m,n);T6(2,1,m,n),T6(2,2,m,n),T6(2,3,m,n);T6(3,1,m,n),T6(3,2,m,n),T6(3,3,m,n)];
        
        T22=[T6(4,4,m,n),T6(4,5,m,n),T6(4,6,m,n);T6(5,4,m,n),T6(5,5,m,n),T6(5,6,m,n);T6(6,4,m,n),T6(6,5,m,n),T6(6,6,m,n)];
        
        V12=[T6(1,4,m,n),T6(1,5,m,n),T6(1,6,m,n);T6(2,4,m,n),T6(2,5,m,n),T6(2,6,m,n);T6(3,4,m,n),T6(3,5,m,n),T6(3,6,m,n)];
        
        dist=zeros(1800,1);%存储距离
        y=zeros(1800,2);
        yy=zeros(1800,2);
        
        %@2 相干区域采样间隔
        for k=1:1800
            th=pi/1800*k;
            T=(T11+T22)/2;
            A=(exp(th*1i)*V12+exp(-th*1i)*V12')/2;
            W=-T\A;
            [d dd]=eig(W);%d是特征向量，dd是特征值
            
            %@2计算Re(y*exp(faio*i))的极大小值，找出对应的特征向量w1,w2,得到相干区域边界点
            re(1)=d(:,1)'*A*d(:,1)/(d(:,1)'*T*d(:,1))*exp(th*1i);
            re(2)=d(:,2)'*A*d(:,2)/(d(:,2)'*T*d(:,2))*exp(th*1i);
            re(3)=d(:,3)'*A*d(:,3)/(d(:,3)'*T*d(:,3))*exp(th*1i);
            re=real(re);

            %判断极大值极小值对应的特征值点，计算两点之间的距离记录下
            re_max=find(re==max(re));
            re_min=find(re==min(re));
            y1=d(:,re_max)'*V12*d(:,re_max)/(d(:,re_max)'*T*d(:,re_max));
            y2=d(:,re_min)'*V12*d(:,re_min)/(d(:,re_min)'*T*d(:,re_min));
             
            dist(k)=abs(y1-y2);
            y(k,1)=angle(y1);y(k,2)=abs(y1);
            yy(k,1)=angle(y2);yy(k,2)=abs(y2);
            
        end
        
        
        one=ones(1,629);
        figure;polar(0:0.01:2*pi,one);hold on;
        polar(y(:,1)',y(:,2)','b.');hold on;
        polar(yy(:,1)',yy(:,2)','y.');hold on;
         
        %@3 计算最大距离对应长轴的两交点，接着计算长轴直线的参数斜距M和截距C，y=Mx+C
        num=find(dist==max(dist));
        thth=pi/1800*num;
        
        T=(T11+T22)/2;
        A=(exp(thth*1i)*V12+exp(-thth*1i)*V12')/2;
        W=-pinv(T)*A;
        [d dd]=eig(W);%d是特征向量，dd是特征值
        
        re(1)=d(:,1)'*A*d(:,1)/(d(:,1)'*T*d(:,1))*exp(th*1i);
        re(2)=d(:,2)'*A*d(:,2)/(d(:,2)'*T*d(:,2))*exp(th*1i);
        re(3)=d(:,3)'*A*d(:,3)/(d(:,3)'*T*d(:,3))*exp(th*1i);
        re=real(re);   

        re_max=find(re==max(re));
        re_min=find(re==min(re));
        Rmax=d(:,re_max)'*V12*d(:,re_max)/(d(:,re_max)'*T*d(:,re_max));
        Rmin=d(:,re_min)'*V12*d(:,re_min)/(d(:,re_min)'*T*d(:,re_min));
        
        polar(angle(Rmax),abs(Rmax),'k+');hold on;
        text(real(Rmax),imag(Rmax),'\color{black} Rmax');
        polar(angle(Rmin),abs(Rmin),'k+');hold on;
        text(real(Rmin),imag(Rmin),'\color{black} Rmin');
        polar(angle(ypdh(m,n)),abs(ypdh(m,n)),'k+');hold on;
        text(real(ypdh(m,n)),imag(ypdh(m,n)),'\color{black} PDHigh');
        polar(angle(ypdl(m,n)),abs(ypdl(m,n)),'k+');hold on;
        text(real(ypdl(m,n)),imag(ypdl(m,n)),'\color{black} PDLow');
        
        b=tan(angle(Rmax-Rmin));
        a=imag(Rmax)-(imag(Rmax)-imag(Rmin))/(real(Rmax)-real(Rmin))*real(Rmax);
        %a=imag(Rmax(m,n))-(imag(Rmax(m,n))-imag(Rmin(m,n)))/(real(Rmax(m,n))-real(Rmin(m,n)))*real(Rmax(m,n));
        
        %@5 绘制长轴直线与地表相位点
        
        nn=1;
        for l=-1:0.001:1
            
            if(l^2+(b*l+a)^2<=1)
                X(nn)=l;
                Y(nn)=b*X(nn)+a;
                nn=nn+1;
            else
                
            end
        end
        
        plot(X,Y,'r');
        
        %@4 计算相干区域长轴与复平面单位圆的交点，并判断地表相位
        %计算拟合直线y=bx+a与单位圆x^2+y^2=1的两个交点
        delta=4*a^2*b^2-4*(1+b^2)*(a^2-1);
        fai1_x=(-2*a*b+sqrt(delta))/(2*(1+b^2));
        fai1_y=b*fai1_x+a;
        fai1=fai1_x+fai1_y*1i;
        fai2_x=(-2*a*b-sqrt(delta))/(2*(1+b^2));
        fai2_y=b*fai2_x+a;
        fai2=fai2_x+fai2_y*1i;
        %将某一交点与HV（体散射占优，代表植被相位）距离大于离HH-VV（二面角散射占优，代表地表相位）的距离比较，为判断地表相位点
        dist_pdh=abs(fai1-ypdh(m,n));
        dist_pdl=abs(fai1-ypdl(m,n));
%           angle1=angle(fai1);
%           angle2=angle(fai2);
        if(dist_pdh>dist_pdl)
            fai0=angle(fai1);
        else
            fai0=angle(fai2);
        end
        polar(angle(exp(fai0*1i)),abs(exp(fai0*1i)),'k+');hold on;
        text(real(exp(fai0*1i)),imag(exp(fai0*1i)),'\color{black} Fai0');
        

        
        
        
    %C 参数输入错误提示
    else
        disp('Parameter format input errors');
    end
end




    