% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File     : function CoherenceSetLoci.m
% Project  : Coherence Set Draw
% Authors  : Gao Han
% Version  : 1.0
% Creation : 2017/9/29
% Update   : 2017/9/29
% *-------------------------------------------------------------------------------
% Description : For drawing the coherence set using classical methods.
%
% Preference:
% (1) Page33 in the Doctoral Paper of Luo HuanMin (Method 1---'Point')
% (2) On The Role of Coherence Optimization in Polarimetric SAR Interferometry,
% by Romeo Tatsambon Fomena (Method 2---'CLM')
% (3) An Interferometric Coherence Optimization Method in Radar Polarimetry for
% High-Resolution Imagery, by Elise Colin (Method 3---'NR')
% Inputs  : You can have 5 variates
% (1)T6: The cell for T6.
% (2)row: The row location of the selected pixel.
% (3)col: The col location of the selected pixel.
% (4)method: The basic method for calculation of the coherence set.('Point'/'CLM'/'NR')
%   ①'Point': Using all possible W, draw all the coherence point;
%   ②'CLM': Constrained Lanrange Multipliers, detailed by the PPT
%   'Description of function CoherenceSetLoci'.
%   ③'NR': Numerical Radius, detailed by the PPT 'Description of
%   function CoherenceSetLoci'.
% (5)isdraw: Whether draw the coherence shape figure or not.(1/0)
%
% Outputs :
% (1)CoPoint: The set of Coherence Set Points by the specific method.
%
% i.e: [Copoint] = CoherenceSetLoci(T6,451,201,'NR',1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CoPoint] = CoherenceSetLoci(T6,row,col,method,isdraw)

i = sqrt(-1);% The complex signal.

if iscell(T6)

% read the T11/T22/T12 of pixel from the cell T6
T11 = [T6{1,1}(row,col),T6{1,2}(row,col),T6{1,3}(row,col);
    T6{2,1}(row,col),T6{2,2}(row,col),T6{2,3}(row,col);
    T6{3,1}(row,col),T6{3,2}(row,col),T6{3,3}(row,col)];

T22 = [T6{4,4}(row,col),T6{4,5}(row,col),T6{4,6}(row,col);
    T6{5,4}(row,col),T6{5,5}(row,col),T6{5,6}(row,col);
    T6{6,4}(row,col),T6{6,5}(row,col),T6{6,6}(row,col)];

T12 = [T6{1,4}(row,col),T6{1,5}(row,col),T6{1,6}(row,col);
    T6{2,4}(row,col),T6{2,5}(row,col),T6{2,6}(row,col);
    T6{3,4}(row,col),T6{3,5}(row,col),T6{3,6}(row,col)];
else
    T11 = [T6(1,1),T6(1,2),T6(1,3);
    T6(2,1),T6(2,2),T6(2,3);
    T6(3,1),T6(3,2),T6(3,3)];

T22 =  [T6(4,4),T6(4,5),T6(4,6);
    T6(5,4),T6(5,5),T6(5,6);
    T6(6,4),T6(6,5),T6(6,6)];

T12 = [T6(1,4),T6(1,5),T6(1,6);
    T6(2,4),T6(2,5),T6(2,6);
    T6(3,4),T6(3,5),T6(3,6)];
end

T =( T11 + T22 ) / 2;% T11 is similar to the T22 when the incidences are the same.
A = sqrtm(pinv(T)) * T12 * sqrtm(pinv(T));

CoPoint = [];% Initilize the coherence points.
m = 200;% 迭代计算的角度步长

switch method
    case 'Point'
        temp = 0.3;
        j=1;
        for alpha = 0:temp:2*pi
            for beita = 0:temp:pi
                for gamma = - pi:temp:pi
                    for yita = -pi:temp:pi;
                        W = [cos(alpha);sin(alpha)*cos(beita)*exp(i*yita);sin(alpha)*sin(beita)*exp(i*gamma)];
                        CoPoint(j) = W' * A * W;
                        j = j + 1;
                    end
                end
            end
        end
        
    case 'CLM'
        for  j = 1 : m+1
            theta = (j-1) * 2 * pi /m;
            TT12 = ( exp(i*theta)*T12 + exp(-i*theta)*T12' ) / 2;
            [V,L]= eig( T\TT12 );
            [Y,Inde] = max(diag(real(L)));
            E = V(:,Inde) / norm(V(:,Inde));
            Tmp12 = E' * T12 * E;
            Tmp11 = E' * T * E;            
            CoPoint(j) = Tmp12 / Tmp11  ;
        end
        
    case 'NR'
        for  j = 1:m+1
            theta = (j-1) * 2 * pi /m;
            Atheta = exp(i*theta)*A;
            H = (Atheta + Atheta')/2;
            [X,L]= eig(H);
            [Y,Inde]  = max(diag(real(L)));
            E = X(:,Inde) / norm(X(:,Inde));
            S = E' * A * E;
            CoPoint(j) = S;
        end
end

if(isdraw)
    switch method 
        case 'Point'
        figure;polar(0,1);hold on
        polar(angle(CoPoint),abs(CoPoint),'+');hold on;
        temp = strcat(method,'CoherenceSet');
        title(temp);
        otherwise
        figure;polar(0,1);hold on
        polar(angle(CoPoint),abs(CoPoint));hold on;
        temp = strcat(method,'CoherenceSet');
        title(temp);
    end
end

% polar(angle(q),abs(q));hold on;
% polar(angle(HH(row,col)),abs(HH(row,col)),'k+');hold on;
% polar(angle(VV(row,col)),abs(VV(row,col)),'b+');hold on;
% polar(angle(HV(row,col)),abs(HV(row,col)),'r+');hold on;
%
% polar(angle(AW1(row,col)),abs(AW1(row,col)),'bo');hold on;
% polar(angle(AW2(row,col)),abs(AW1(row,col)),'r*');hold on;
% polar(angle(AW3(row,col)),abs(AW3(row,col)),'g+');hold on;
