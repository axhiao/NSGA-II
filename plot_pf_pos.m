position = 'E:\进化多目标方法工具箱\NSGAII-DE-Store';

str_pf = strcat(position,'\POF\');
str_pos = strcat(position,'\POS\');
for function_num = 1:1  
    for c=1:50
        %
        data_pf = strcat(str_pf,num2str(c));
        data_pf = strcat(data_pf,num2str(0));
        data_pf = strcat(data_pf,'\POF_NSGA2_UF');
        data_pf = strcat(data_pf,num2str(function_num));
        data_pf = strcat(data_pf,'_RUN1.txt');
        pf = importdata(data_pf);
        %
        data_pos = strcat(str_pos,num2str(c));
        data_pos = strcat(data_pos,num2str(0));
        data_pos = strcat(data_pos,'\POS_NSGA2_UF');
        data_pos = strcat(data_pos,num2str(function_num));
        data_pos = strcat(data_pos,'_RUN1.txt');
        pos = importdata(data_pos);    
    %% show
        switch function_num
            case {1 2 3}
                x=(0:0.01:1)';
                y=1-x.^0.5;
                subplot(1,2,1);
                plot(pf(:,1),pf(:,2),'ro',x(:),y(:),'-g*');
                axis([0 1.1 0 1.1]);
                
                subplot(1,2,2);
                switch function_num
                    case 1
                        x=(0:0.001:1)';
                        xx=ones(length(x),3);
                        xx(:,1)=x;
                        for i=1:length(x)
                            for j=2:3
                                xx(i,j)=sin(6*pi*xx(i,1)+j*pi/30);
                            end
                        end
                        plot3(pos(:,1),pos(:,2),pos(:,3),'go',xx(:,1),xx(:,2),xx(:,3),'-b');
                        xlabel('x1');
                        ylabel('x2');
                        zlabel('x3');
                        title('UF1-PS');
                    case 2
                        x=(0:0.001:1)';
                        xx=ones(length(x),3);
                        xx(:,1)=x;
                        for i=1:length(x)
                            for j=2:3
                                if mod(j,2)==1
                                    xx(i,j)=(0.3*xx(i,1)*xx(i,1)*cos(24*pi*xx(i,1)+4*j*pi/30)+0.6*xx(i,1))*cos(6*pi*xx(i,1)+j*pi/30);
                                else
                                    xx(i,j)=(0.3*xx(i,1)*xx(i,1)*cos(24*pi*xx(i,1)+4*j*pi/30)+0.6*xx(i,1))*sin(6*pi*xx(i,1)+j*pi/30);
                                end
                            end
                        end
                        plot3(pos(:,1),pos(:,2),pos(:,3),'go',xx(:,1),xx(:,2),xx(:,3),'-b');
                        xlabel('x1');
                        ylabel('x2');
                        zlabel('x3');
                        title('UF2-PS');
                    case 3
                        x=(0:0.001:1)';
                        xx=ones(length(x),3);
                        xx(:,1)=x;
                        for i=1:length(x)
                            for j=2:3
                                xx(i,j)=xx(i,1)^(0.5*(1+3.0*(j-2)/(30-2)));
                            end
                        end
                        plot3(pos(:,1),pos(:,2),pos(:,3),'go',xx(:,1),xx(:,2),xx(:,3),'-b');
                        xlabel('x1');
                        ylabel('x2');
                        zlabel('x3');
                        title('UF3-PS');
                    otherwise
                end
            case 4
                x=(0:0.01:1)';
                y=1-x.^2;
                subplot(1,2,1);
                plot(pf(:,1),pf(:,2),'ro',x(:),y(:),'-g*')
                xlabel('f1');
                ylabel('f2');
                title('UF4');
                h=legend('MOEA-D-DE result','真实的Pareto前端',2);

                subplot(1,2,2);
                x=(0:0.001:1)';
                xx=ones(length(x),3);
                xx(:,1)=x;
                for i=1:length(x)
                    for j=2:3
                        xx(i,j)=sin(6*pi*xx(i,1)+j*pi/30);
                    end
                end
                plot3(pos(:,1),pos(:,2),pos(:,3),'go',xx(:,1),xx(:,2),xx(:,3),'-b')
                xlabel('x1');
                ylabel('x2');
                zlabel('x3');
                title('UF4-PS');
            case 5
                N = 10;
                x=(0:1/(2*N):1);
                y=1-x;
                subplot(1,1,1);
                plot(pf(:,1),pf(:,2),'ro',x(:),y(:),'g*')
                xlabel('f1');
                ylabel('f2');
                title('UF5');
                h=legend('MOEA-D-DE result','真实的Pareto前端',2);
            case 6
                N = 2;
                x1=(1/4:0.005:0.5);
                y1=1-x1;
                x2=(3/4:0.005:1);
                y2=1-x2;
                subplot(1,1,1);
                plot(pf(:,1),pf(:,2),'ro',x1(:),y1(:),'-g*',x2(:),y2(:),'-g*',0,1,'g*')
                xlabel('f1');
                ylabel('f2');
                title('UF6');
                h=legend('MOEA-D-DE result','真实的Pareto前端',2);
            case 7
                x=(0:0.01:1)';
                y=1-x;
                subplot(1,2,1);
                plot(pf(:,1),pf(:,2),'ro',x(:),y(:),'-g*')
                xlabel('f1');
                ylabel('f2');
                title('UF7');
                h=legend('MOEA-D-DE result','真实的Pareto前端',2);
                
                x=(0:0.001:1)';
                xx=ones(length(x),3);
                xx(:,1)=x;
                for i=1:length(x)
                    for j=2:3
                        xx(i,j)=sin(6*pi*xx(i,1)+j*pi/30);
                    end
                end
                subplot(1,2,2);
                plot3(pos(:,1),pos(:,2),pos(:,3),'go',xx(:,1),xx(:,2),xx(:,3),'-b')
                xlabel('x1');
                ylabel('x2');
                zlabel('x3');
                title('UF7-PS');
            otherwise
        end
        pause(1);
    end
%     data_pf = strcat(str_pf,'\POF_NSGA2_UF');
%     data_pf = strcat(data_pf,num2str(function_num));
%     data_pf = strcat(data_pf,'_RUN1.txt');
%     pf = importdata(data_pf);
%     
%     data_pos = strcat(str_pos,'\POS_NSGA2_UF');
%     data_pos = strcat(data_pos,num2str(function_num));
%     data_pos = strcat(data_pos,'_RUN1.txt');
%     pos = importdata(data_pos);    
end