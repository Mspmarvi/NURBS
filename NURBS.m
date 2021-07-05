classdef NURBS<handle
    %Input Argument Description
    %1- at first you must build a matrix that it's size must be (n+1,3).it
    %means we have two vectors integrated in one matrix, first vector is
    %for x coordinate and second vector is for y coordinate an for z coordinates,this matrix is
    %called controlpoint.
    %2- you can add a number for degree
    %3- at third, we must enter loop,for this part we must add 'close' or
    %'open',nothing more
    %4- at fourth, if you want a B-Spline curve you must add number 1
    %unless you must add vector of weight  with size of (n+1,1)
    %5- at fifth,(it's about boundry condition like uniform or clamped)
    %if you don't want to add anything, the curve have uniform boundry
    %condition else you must add clamped or uniform.
    %   Detailed explanation goes here
    properties
        Number_of_ControlPoints
        Degree
        ControlPoints
        Boundry
        Loop
        Weights
        UVector
        CurveData
        DataPoint
    end   
    methods
        function obj=NURBS(controlpoint,degree,loop,weight,varargin) % this is the constructor of NURBS curve;
%controlpoint is a vector with size (n+1,3) in which n+1 is number of control points ;
% degree is scalar; loop : 'open' or 'close';
% weight is a vector with size (n+1,1), !(if you want to build a B-spline curve just add 1); 
%varargin: if you want to have uniform knot vector don't fill it but if you want clamped curve just add 'clamped' 
            %check inputs
            if nargin<4
                error('not enough inputs')
            elseif nargin>5
                error('too many inputs')
            elseif nargin==4
                boundry='uniform'  ;  %default boundry
            else
                boundry=varargin{1};
            end
            
            %degree
            if isreal(degree)
                obj.Degree=degree;
            else
                error('wrong degree')
            end
            
            %loop
            if ischar(loop)
                obj.Loop=loop;
            else
                error('wrong loop')
            end
            
            %to find numbers of controlpoints
            n=size(controlpoint,1)-1; %number of controlpoints
            obj.Number_of_ControlPoints=n+1;
                        
            if isreal(controlpoint)
                obj.ControlPoints=controlpoint;
            else
                error('wrong controlpoints');
            end
            
            
            if ~isreal(weight)
                error('wrong weight')
            end
            
            %loop
            switch loop
                case 'open'
                    
                    if ischar(boundry)
                        obj.Boundry=boundry;
                    else
                        error('wrong boundry condition')
                    end
                    
                case 'close'
                    CP=controlpoint;
                    controlpoint=zeros(n+1+degree,3);
                    controlpoint(1:n+1,:)=CP(:,:);
                    controlpoint(n+1+1:end,:)=CP(1:degree,:);
                    boundry='uniform';
                    obj.Boundry=boundry;
                    
                otherwise
                    error('you were wrong in putting loop')
            end
            
            n=size(controlpoint,1)-1;
            m=n+degree+1;
            U=zeros(m+1,1);  %list of knot vectors
        
            switch boundry
                case 'uniform'                %we assumed that curve is uniform
                    for i=1:m
                        U(i+1)=U(i)+1/m;
                    end
                case 'clamped'
                    for i=1+degree:m-degree
                        U(i+1)=U(i)+1/(m-2*degree);
                    end
                    U(m+1-degree:m+1)=1;
            end
            obj.UVector=U;
            %to impose the weights
            
            if weight==1
                weight=ones(n+1,1);
                obj.Weights=weight;
            else
                obj.Weights=weight;
            end
            
            if strcmp(loop,'close')
                W=zeros(n+1,1);
                W(1:size(weight,1),1)=weight(:);
                W(size(weight,1)+1:end,1)=weight(1:obj.Degree,1);
                weight=W;
            elseif size(weight,1)~=size(controlpoint,1)
                error('wrong weight or control point, you should check the size of them')
            end
            %%
            obj = Cal_CurveData(obj);
            %%
        end
        function [x,y,z]=Calculate(obj,u_parameter) % to caluclate coordiante of NURBS curve in parameter  u_parameter;  0<= u_parameter<=1 
            x=0; y=0; z=0;w=0;
            B=NURBS.BasisFunction(u_parameter,obj.Degree,obj.UVector);
            for i=1:obj.Number_of_ControlPoints
                x=x+B(i,end)*obj.ControlPoints(i,1);
                y=y+B(i,end)*obj.ControlPoints(i,2);
                z=z+B(i,end)*obj.ControlPoints(i,3);
                w=w+B(i,end)*obj.Weights(i);  
            end
            x=x/w;
            y=y/w;
            z=z/w;
        end
        function plot(obj)
            obj = Cal_CurveData(obj);            
            %             figure(1)
            hold on
            plot3(obj.CurveData(:,1),obj.CurveData(:,2),obj.CurveData(:,3),'b');
            %%
            breakpoint=[];
            m = obj.Number_of_ControlPoints+obj.Degree;
            for i=1+obj.Degree:m-obj.Degree
                u=obj.UVector(i);
                [x,y,z]=Calculate(obj,u);
                breakpoint(end+1,1)=x;
                breakpoint(end,2)=y;
                breakpoint(end,3)=z;
            end
            
            %% plotting
            axis equal
            plot3(breakpoint(:,1),breakpoint(:,2),breakpoint(:,3),'mo','MarkerFaceColor',[1 0 0.5],'markersize',1)
            plot3(obj.CurveData(end,1),obj.CurveData(end,2),obj.CurveData(end,3),'mo','MarkerFaceColor',[1 0 0.5],'markersize',1)
            plot3(obj.CurveData(:,1),obj.CurveData(:,2),obj.CurveData(:,3),'r','linewidth',1)
            plot3(obj.ControlPoints(:,1),obj.ControlPoints(:,2),obj.ControlPoints(:,3),'b','linewidth',1);
            plot3(obj.ControlPoints(:,1),obj.ControlPoints(:,2),obj.ControlPoints(:,3),'mo','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'markersize',2);
            switch obj.Loop
                case 'open'
                    for i=1:obj.Number_of_ControlPoints
                        %                          text(controlpoint(i,1),controlpoint(i,2),controlpoint(i,3),['P',num2str(i-1)],'fontsize',9,'HorizontalAlignment','center')
                    end
                case 'close'
                    for i=1:obj.Number_of_ControlPoints
                        %                         text(controlpoint(i,1),controlpoint(i,2),controlpoint(i,3),num2str(i-1),'fontsize',13)
                    end
            end
            
        end
        function plotBasisFunction(obj)
            figure(2)
            hold on
            for i=1:obj.Number_of_ControlPoints
                basis=[];
                for u=linspace(obj.UVector(obj.Degree+1),obj.UVector(end-obj.Degree),700)
                    B=BasisFunction(u,obj.Degree,obj.UVector);
                    basis(end+1)=B(i,end);
                end
                plot(linspace(obj.UVector(obj.Degree+1),obj.UVector(end-obj.Degree),700),basis,'b');
            end
            axis equal
            axis([0 1 0 1])
        end
        function KnotInsertion(obj,uinsert)% if you want to insert a knot with shape preserving; 0<=uinsert<=1
            
            mplus=obj.Number_of_ControlPoints+obj.Degree+1;
            %searching the location of uinsert
            for i=1:mplus-1
                if uinsert>=obj.UVector(i) && uinsert<obj.UVector(i+1)
                    k=i;
                    break
                end
            end
            
            %calculate a(i) you must check your cadcam course
            %remember the index here is different from the cadcam course
            
            a=zeros(obj.Number_of_ControlPoints+1,1);
            
            for i=k-obj.Degree+1:k
                a(i)=(uinsert-obj.UVector(i))/(obj.UVector(i+obj.Degree)-obj.UVector(i));
            end
            
            %homogenous space
            H=zeros(obj.Number_of_ControlPoints,3);
            
            for i=1:obj.Number_of_ControlPoints
                H(i,1)=obj.ControlPoints(i,1)*obj.Weights(i,1);
                H(i,2)=obj.ControlPoints(i,2)*obj.Weights(i,1);
                H(i,3)=obj.ControlPoints(i,3)*obj.Weights(i,1);
            end
            %to calculate the new homogenous control points & weight
            NH=zeros(obj.Number_of_ControlPoints+1,3);
            NW=zeros(obj.Number_of_ControlPoints+1,1);
            NCP=zeros(obj.Number_of_ControlPoints+1,3);
            
            for i=k-obj.Degree+1:k
                NH(i,1)=(1-a(i))*H(i-1,1)+a(i)*H(i,1);
                NH(i,2)=(1-a(i))*H(i-1,2)+a(i)*H(i,2);
                NH(i,3)=(1-a(i))*H(i-1,3)+a(i)*H(i,3);
                
                NW(i,1)=(1-a(i))*obj.Weights(i-1,1)+a(i)*obj.Weights(i,1);
                
                NCP(i,1)=NH(i,1)/NW(i,1);
                NCP(i,2)=NH(i,2)/NW(i,1);
                NCP(i,3)=NH(i,3)/NW(i,1);
            end
            NCP(1:k-obj.Degree,:)=obj.ControlPoints(1:k-obj.Degree,:);
            NCP(k+1:end,:)=obj.ControlPoints(k:end,:);
            obj.ControlPoints=NCP;
            %%%
            NW(1:k-obj.Degree,:)=obj.Weights(1:k-obj.Degree,:);
            NW(k+1:end,:)=obj.Weights(k:end,:);
            obj.Weights=NW;
            %%%
            DU=zeros(mplus+1,1);
            DU(1:k,1)=obj.UVector(1:k,1);
            DU(k+1,1)=uinsert;
            DU(k+2:end,1)=obj.UVector(k+1:end,1);
            obj.UVector=DU;
            %%%
            obj.Number_of_ControlPoints=size(obj.ControlPoints,1);
        end
        function obj = Interpolation(obj) %to use interpolation you can use datapoint as controlpoint to build the nurbs curve then use this function to interpolate controlpoint as datapoint
                
                X=obj.ControlPoints(:,1);
                Y=obj.ControlPoints(:,2);
                Z=obj.ControlPoints(:,3);
                
                obj.DataPoint=[X,Y,Z];
                
                sizedata=size(X);
                nplus=sizedata(1,1);
                n=nplus-1;
                p=obj.Degree;
                m=p+n+1;
                t=zeros(n+1,1);
                
                %this method is chord length
%                 L=0;
%                 for i=1:n
%                     L=L+sqrt((X(i+1)-X(i))^2+(Y(i+1)-Y(i))^2+(Z(i+1)-Z(i))^2);
%                 end
%                 for i=1:n
%                     t(i+1)=t(i)+sqrt((X(i+1)-X(i))^2+(Y(i+1)-Y(i))^2+(Z(i+1)-Z(i))^2);
%                 end
%                 t=t/L;
                %%this is uniform distribution for t(k)
                for i=1:n
                    t(i+1)=i/n;
                end
                
                U=zeros(m+1,1);
                U(m-p+1:m+1)=1;
                for i=1:m-2*(p+1)+1
                    U(i+p+1)=sum(t(i+1:i+p))/p;
                end
                NN=zeros(n+1,n+1);
                N=zeros(m,p+1);
                for i=1+p:m-p
                    N(:,1)=0;
                    N(i,1)=1;
                    for d=1:n+1
                        if t(d)<U(i+1) && t(d)>=U(i)
                            for q=1:p
                                for j=1:m-q
                                    left=(t(d)-U(j))/(U(j+q)-U(j));
                                    if isinf(left)==1 || isnan(left)==1
                                        left=0;
                                    end
                                    right=(U(j+q+1)-t(d))/(U(j+q+1)-U(j+1));
                                    if isinf(right)==1 || isnan(right)==1
                                        right=0;
                                    end
                                    N(j,q+1)=left.*N(j,q)+right.*N(j+1,q);
                                end
                            end
                            for k=1:n+1
                                NN(d,k)=N(k,p+1);
                            end
                        end
                    end
                end
                
                NN(n+1,n+1)=1;
                Px=inv(NN)*X;
                Py=inv(NN)*Y;
                Pz=inv(NN)*Z;
                                                
                obj.ControlPoints=[Px,Py,Pz];
                obj.Number_of_ControlPoints=size(obj.ControlPoints,1);
                obj.UVector=U;
                obj.Boundry='clamped';
                obj.Loop='open';
                obj = Cal_CurveData(obj);
        end
        function [dx,dy,dz] = B_spline_derivative(obj,u_parameter,nth_derivative) % to calculate derivation nth of NURBS 
            
            CP=obj.ControlPoints;
            n=obj.Number_of_ControlPoints-1;
            Q=CP;
            nth=nth_derivative;
            Uvector=obj.UVector;
            for i=1:nth
                newQ=zeros(n+1-i,3);
                for j=1:n+1-i
                    newQ(j,1)=(obj.Degree-i+1)*(Q(j+1,1)-Q(j,1))/(Uvector(j+obj.Degree+1)-Uvector(j+i));
                    newQ(j,2)=(obj.Degree-i+1)*(Q(j+1,2)-Q(j,2))/(Uvector(j+obj.Degree+1)-Uvector(j+i));
                    newQ(j,3)=(obj.Degree-i+1)*(Q(j+1,3)-Q(j,3))/(Uvector(j+obj.Degree+1)-Uvector(j+i));
                end
                Q=newQ;
            end
            Uvector=Uvector(nth+1:end-nth);
            Basis=BasisFunction(u_parameter,obj.Degree,Uvector);
            
            dx=0; dy=0;dz=0;
            for i=1:size(Q,1)
                dx=dx+Basis(i,end-nth)*Q(i,1);
                dy=dy+Basis(i,end-nth)*Q(i,2);
                dz=dz+Basis(i,end-nth)*Q(i,3);
            end
            
        end
        function plotDerivates(obj)
            figure(1)
            hold on
            for u=0:0.01:1
                [x,y,z]=Calculate(obj,u);
                [dx,dy,dz]=B_spline_derivative(obj,u,1);
                dp=[dx,dy,dz];
                dp=dp/10;
                P=[x,y,z];
                dP=[x+dp(1),y+dp(2),z+dp(3)];
                plot3([P(1),dP(1)],[P(2),dP(2)],[P(3),dP(3)],'k');
            end
        end
        function obj = Approximation(obj,distribution,division,varargin)% same as interpolation function; varargin is number of controlpoints if you want to add
          m=obj.Number_of_ControlPoints-1;
          obj.DataPoint=obj.ControlPoints;
          D=obj.ControlPoints;
          p=obj.Degree;
          tk=zeros(m+1,1);
          switch distribution
              case  'uniform'
                  for i=1:m
                      tk(i+1)=tk(i)+1/(m);
                  end
              case   'chord'
                  L=0;
                  for i=1:m
                      L=L+sqrt((D(i+1,1)-D(i,1))^2+(D(i+1,2)-D(i,2))^2+(D(i+1,3)-D(i,3))^2);
                  end
                  for i=1:m
                      tk(i+1)=tk(i)+sqrt((D(i+1,1)-D(i,1))^2+(D(i+1,2)-D(i,2))^2+(D(i+1,3)-D(i,3))^2)/L;
                  end
              case    'centripetal'
                  L=0;
                  for i=1:m
                      L=L+(sqrt((D(i+1,1)-D(i,1))^2+(D(i+1,2)-D(i,2))^2+(D(i+1,3)-D(i,3))^2))^0.5;
                  end
                  for i=1:m
                      tk(i+1)=tk(i)+((sqrt((D(i+1,1)-D(i,1))^2+(D(i+1,2)-D(i,2))^2+(D(i+1,3)-D(i,3))^2))^0.5)/L;
                  end
              otherwise
                  error('wrong statement')
          end
          %%
          if numel(varargin)==0
              n=floor((m+1)/division)+1;
          else
              n=varargin{1}-1;
          end
          mm=n+p+1;
          du=1/mm;
          U=zeros(mm+1,1);
          for i=1:mm
              U(i+1)=U(i)+du;
          end
          U(1:p+1)=0;
          U(mm-p+1:mm+1)=1;
           %%
           N=zeros(m-1,n-1);
           for i=1:m-1
               B=BasisFunction(tk(i+1),p,U);
               for j=1:n-1
                   N(i,j)=B(j+1,end);
               end
           end
           %%
           Rk=zeros(m-1,3);
           for i=1:size(Rk,1)
               B=BasisFunction(tk(i+1),p,U);
               Rk(i,1)=D(i+1,1)-B(1,end)*D(1,1)-B(n+1,end)*D(end,1);
               Rk(i,2)=D(i+1,2)-B(1,end)*D(1,2)-B(n+1,end)*D(end,2);
               Rk(i,3)=D(i+1,3)-B(1,end)*D(1,3)-B(n+1,end)*D(end,3);
           end
           %%
           R=zeros(n-1,3);
           for i=1:n-1
               rx=0; ry=0; rz=0;
              for j=1:m-1
                  B=BasisFunction(tk(j+1),p,U);
                  rx=rx+B(i+1,end)*Rk(j,1);
                  ry=ry+B(i+1,end)*Rk(j,2);
                  rz=rz+B(i+1,end)*Rk(j,3);
              end
              R(i,1)=rx; R(i,2)=ry; R(i,3)=rz;
           end
           P=inv(N'*N)*R;
           PP=zeros(size(P,1)+2,3);
           PP(1,:)=D(1,:);
           PP(end,:)=D(end,:);
           PP(2:end-1,:)=P(:,:);
           obj.ControlPoints=PP;
           obj.UVector=U;
           
           obj.Number_of_ControlPoints=n+1;
           obj.Weights = ones(n+1,1);
           obj = Cal_CurveData(obj);
        end
        function plotDataPoint(obj)
            plot3(obj.DataPoint(:,1),obj.DataPoint(:,2),obj.DataPoint(:,3),'.','markerfacecolor','k','markeredgecolor','k')
        end
        function obj = Close(obj) % convert open curve to close curve
            newcp = obj.ControlPoints(1:obj.Degree+1,:);
            newweight = obj.Weights(1:obj.Degree+1,1);
            cp = [obj.ControlPoints
                newcp];
            w = [obj.Weights
                newweight];
            
            obj.ControlPoints = cp; obj.Weights = w; obj.Loop = 'close'; obj.Number_of_ControlPoints = size(obj.ControlPoints,1);
            
            m=obj.Degree+obj.Number_of_ControlPoints;
            obj.UVector=linspace(0,1,m+1);
        end
        function obj = clamped2uniform(obj)% convert clamped curve into uniform curve
            n = obj.Number_of_ControlPoints-1;
            m=n+obj.Degree+1;
            U=zeros(m+1,1);  %list of knot vectors
               %we assumed that curve is uniform
                    for i=1:m
                        U(i+1)=U(i)+1/m;
                    end
            obj.UVector=U;
            obj = Cal_CurveData(obj);
        end
        function obj = Cal_CurveData(obj) % calculation of whole curve 
            number_of_points = 500;
            if strcmp(obj.Loop,'close')
                u = linspace(obj.UVector(obj.Degree+1),obj.UVector(end-obj.Degree-(1)),number_of_points);
                points = zeros(numel(u),3);
                for i = 1:numel(u)
                    [x,y,z] = Calculate(obj,u(i));
                    points(i,:) = [x,y,z];
                end
                obj.CurveData = points;
            else
                u = linspace(obj.UVector(obj.Degree+1),obj.UVector(end-obj.Degree),number_of_points);
                points = zeros(numel(u),3);
                for i = 1:numel(u)
                    [x,y,z] = Calculate(obj,u(i));
                    points(i,:) = [x,y,z];
                end
                obj.CurveData = points;
            end
        end
        function plotCurveData(obj)
            obj = Cal_CurveData(obj);
            plot3(obj.CurveData(:,1),obj.CurveData(:,2),obj.CurveData(:,3),'.k');
        end
    end
    methods(Static)
        function  BasisMatrix=BasisFunction(u_parameter,degree,Uvector)
            m=numel(Uvector)-1;
            p=degree;
            U=Uvector;
            BasisMatrix=zeros(m,p+1);
            u=u_parameter;
            for i=1:m
                if u>=U(i) && u<U(i+1)
                    location=i;
                    break
                end
            end
            if u==1
                location=m-degree;
            end
            BasisMatrix(:,1)=0;
            BasisMatrix(location,1)=1;
            for q=1:p
                for j=1:m-q
                    left=((u-U(j))/(U(j+q)-U(j)));
                    if  isnan(left)==1 || isinf(left)==1
                        left=0;
                    end
                    right=((U(j+q+1)-u)/(U(j+q+1)-U(j+1)));
                    if isnan(right)==1 || isinf(right)==1
                        right=0;
                    end
                    BasisMatrix(j,q+1)=left.*BasisMatrix(j,q)+right.*BasisMatrix(j+1,q);
                end
            end
        end
    end
end