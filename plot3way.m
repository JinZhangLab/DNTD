function plot3way(A,B,C)
% Plot profiles of tree-way data decompostion
% if shiting occurs in B model, the function is also useful.
% A, B, and C are the resovled profiles.
% Written by J zhang
% zhangjin@mail.nankai.edu.cn
DimA=size(A);
DimB=size(B);
DimC=size(C);

[axis_A,axis_B,axis_C]=deal(1:DimA(1),1:DimB(1),1:DimC(1));

figure;
isshift=length(DimB)==3;
for i=1:DimA(2)
    subplot(DimA(2),3,(i-1)*3+1);
    plot(axis_A,A(:,i),'-o','LineWidth',2);xlabel('A');ylabel('Intensity');
    
    subplot(DimA(2),3,(i-1)*3+2);
    if ~isshift
        plot(axis_B,B(:,i),'-','LineWidth',1);xlabel('B');ylabel('Intensity');
    else
        hold on;
        for j=1:DimB(3)
            Bj=squeeze(B(:,:,j));
            plot(axis_B,Bj(:,i),'-','LineWidth',1);xlabel('B');ylabel('Intensity');
        end
    end
    
    subplot(DimA(2),3,(i-1)*3+3);
    plot(axis_C,C(:,i),'-s','LineWidth',2);xlabel('C');ylabel('Intensity');
end
end


