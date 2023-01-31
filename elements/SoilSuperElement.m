classdef SoilSuperElement < Element
    %SOILSUPERELEMENT Base class for soil superelement
    %   Detailed explanation goes here
    
    properties (Access = private)         
        lengthX
        lengthY
        
    end    
    
    methods
        % Constructor
        function soilsuperelement = SoilSuperElement(id, nodeArray, Nx)
            requiredPropertyNames = cellstr(["THICKNESS", "YOUNGS_MODULUS", "POISSON_RATIO", "DENSITY", "ELEMENTAL_DAMPING", ...
                "YOUNGS_MODULUS_HS", "POISSON_RATIO_HS", "DENSITY_HS", "ELEMENTAL_DAMPING_HS", ...
                "SAMPLE_NUMBER", "SAMPLE_RATE", "FREQUENCY"]);
            
             % define the arguments for the super class constructor call
            if nargin == 0
                super_args = {};
            elseif nargin == 3
                if ~(length(nodeArray) == Nx^2 && isa(nodeArray,'Node'))
                    error('problem with the nodes in element %d', id);
                end
                super_args = {id, nodeArray, requiredPropertyNames};
            end
            
            % call the super class constructor
            soilsuperelement@Element(super_args{:});
            soilsuperelement.dofNames = cellstr(["DISPLACEMENT_X", "DISPLACEMENT_Y", ...
                "DISPLACEMENT_Z"]);
             
        end
    
        function barycenter(obj)
        end
        
        %Initialization
        function initialize(obj)
            obj.lengthX = computeLength(obj.nodeArray(1).getCoords, ...
                obj.nodeArray(end).getCoords);
            
            obj.lengthY = obj.lengthX / sqrt(2);
            obj.lengthX = obj.lengthY; 
        end                
        
        function stiffnessmatrix = computeLocalStiffnessMatrix(obj)
            % compute stiffness matrix
            % estimation of material properties for calculation of stiffness matrix
            
            Emodul = obj.getPropertyValue('YOUNGS_MODULUS');
            D = obj.getPropertyValue('ELEMENTAL_DAMPING');
            Emodul = Emodul*(1-2*D*1i);
            
            PoissonRatio = obj.getPropertyValue('POISSON_RATIO');
            rho = obj.getPropertyValue('DENSITY');
            
            h = obj.getPropertyValue('THICKNESS');
            
            G=Emodul./(2*(1+PoissonRatio));
            mue=G;
            lambda=2*G.*PoissonRatio./(1-2*PoissonRatio);
            cp=sqrt((lambda+2*mue)./rho);
            cs=sqrt(mue./rho);
            
            Emodul_HS = obj.getPropertyValue('YOUNGS_MODULUS_HS');
            D_HS = obj.getPropertyValue('ELEMENTAL_DAMPING_HS');
            Emodul_HS = Emodul_HS*(1-2*D_HS*1i);
            
            PoissonRatio_HS = obj.getPropertyValue('POISSON_RATIO_HS');
            rho_HS = obj.getPropertyValue('DENSITY_HS');
            
            G_HS=Emodul_HS./(2*(1+PoissonRatio_HS));
            mue_HS=G_HS;
            lambda_HS=2*G_HS.*PoissonRatio_HS./(1-2*PoissonRatio_HS);
            cp_HS=sqrt((lambda_HS+2*mue_HS)./rho_HS);
            cs_HS=sqrt(mue_HS./rho_HS);
            
            ITM_PROPERTIES = [mue, lambda, cp, cs, h, mue_HS, lambda_HS, cp_HS, cs_HS];
            
            % In order to gain stiffnessMatrix using symmetry characteristics of the halfspace,
            % displacemts have to be evaluated over twice the number of samples
            Nx = 2 * obj.getPropertyValue('SAMPLE_NUMBER');
            dx = obj.getPropertyValue('SAMPLE_RATE');
            Nt = 16;
            frequency = obj.getPropertyValue('FREQUENCY');
            
            Bx=dx*Nx;
            kx0=2*pi/Bx;
            kx=-Nx/2*kx0:kx0:(Nx/2-1)*kx0;
            
            stiffnessmatrix = obj.computeStiffnessMatrixSoil( frequency,Nx,Nt,kx,dx,ITM_PROPERTIES );
        end
        
        function massmatrix = computeLocalMassMatrix(obj)
            Nx = obj.getPropertyValue('SAMPLE_NUMBER');
            massmatrix = sparse (3*(Nx^2),3*(Nx^2));
        end
        
        function dampingmatrix = computeLocalDampingMatrix(obj)
            Nx = obj.getPropertyValue('SAMPLE_NUMBER');
            dampingmatrix = sparse (3*(Nx^2),3*(Nx^2));
        end
        
        function dofs = getDofList(obj)
            len = length(obj.nodeArray)*3;
            
            dofs(1:3:len) = obj.nodeArray.getDof('DISPLACEMENT_X');
            dofs(2:3:len) = obj.nodeArray.getDof('DISPLACEMENT_Y');
            dofs(3:3:len) = obj.nodeArray.getDof('DISPLACEMENT_Z');
            
        end
        
        function vals = getValuesVector(obj, step)
            len = length(obj.nodeArray)*3;
            vals = zeros(1,len);
            
            vals(1:3:len) = obj.nodeArray.getDofValue('DISPLACEMENT_X',step);
            vals(2:3:len) = obj.nodeArray.getDofValue('DISPLACEMENT_Y',step);
            vals(3:3:len) = obj.nodeArray.getDofValue('DISPLACEMENT_Z',step);
        end
        
        function update(obj)
        end
        
        function F = computeLocalForceVector(obj)
            Nx = obj.getPropertyValue('SAMPLE_NUMBER');
            F = zeros(3*(Nx)^2,1);
        end
        
        function nc = getNodeConnectivity(obj)
            %nc is an n x 5 array with [sub element id node_id_1 2 3 4]
            
            Nx = obj.getPropertyValue('SAMPLE_NUMBER');
            dx = obj.getPropertyValue('SAMPLE_RATE');
            
            nc =zeros((Nx-1)^2,5);
            nc(:,1)=1:1:size(nc,1);
            
            x=0:dx:(Nx-1)*dx;
            y=x;
            
            j=0;
            for ny = 1:Nx-1
                for nx = 1:Nx-1
                    j=j+1;
                    nc(j,2)= getLocalNodeId(obj.nodeArray,x(nx),y(ny));
                    nc(j,3)= nc(j,2)+1;
                    nc(j,4)= nc(j,3)+Nx;
                    nc(j,5)= nc(j,4)-1;
                end
            end
            
        end
        
    end
    
    methods (Access = private)
        
        function K = computeStiffnessMatrixSoil( obj, frequency,Nx,Nt,kx,dx,ITM_PROPERTIES )
            % compute stiffness matrix using matrix of displacements U in frequency domain (x,y,z,omega)
            % flexibility matrix F: in each column of F there are displacements resulting from
            % point load at each nodes
            F=zeros(3*((Nx/2)^2));
            
            U = obj.displacementmatrix(frequency,Nx,Nt,kx,dx,ITM_PROPERTIES);
            
            for j=1:Nx/2
                for i=1:Nx/2
                    U_MATRIX=U(3*(Nx/2+1-j)+1:3*(Nx+1-j),3*(Nx/2+1-i)+1:3*(Nx+1-i));
                    
                    F(:,3*((i-1)+(j*Nx/2-Nx/2))+1:3*((i-1)+(j*Nx/2-Nx/2))+3)= ...
                        horzcat(reshape(U_MATRIX(1:3:3*Nx/2,:).',3*(Nx/2)^2,1), ...
                        reshape(U_MATRIX(2:3:3*Nx/2,:).',3*(Nx/2)^2,1), reshape(U_MATRIX(3:3:3*Nx/2,:).',3*(Nx/2)^2,1));
                end
            end
            
            K = eye(size(F))/(F);
        end
        
        
        function U = displacementmatrix(obj,frequency,Nx,Nt,kx,dx,ITM_PROPERTIES)
            % compute flexibility matrix F
            % in each column of F there are displacements resulting from
            % point load at each nodes
            U=zeros(3*Nx);
            
            ux = obj.displacement(frequency,Nx,Nt,kx,dx,ITM_PROPERTIES);
            
            for nx=1:3:3*Nx
                for ny=1:Nx
                    U(nx,(3*ny-2):(3*ny))=ux((nx+2)/3,ny,Nt/2,1:3);
                    U(nx+1,(3*ny-2):(3*ny))=ux((nx+2)/3,ny,Nt/2,4:6);
                    U(nx+2,(3*ny-2):(3*ny))=ux((nx+2)/3,ny,Nt/2,7:9);
                end
            end
        end
        
        
        function ux = displacement(obj,frequency,Nx,Nt,kx,dx,ITM_PROPERTIES)
            % compute displacements at each nodes in frequency domain (x,y,z,omega)
            u=zeros(Nx,Nx,Nt/2+1,9);
            
            ky = kx;
            
            Pt = obj.computeLoadSoilSuperElement(Nx,Nt,dx,frequency);
            
            omega0=2*pi*frequency;
            omega=-Nt/2*omega0:omega0:0;
            
            for nt=1:Nt/2+1
                for ny=1:Nx
                    for nx=1:Nx
                        if omega(nt)==0 && kx(nx)==0 && ky(ny)==0
                            
                            u(nx,ny,nt,:)=obj.displacementsHR_stat(kx(nx-1),ky(ny-1),Pt(nx,ny,nt),ITM_PROPERTIES);
                        elseif omega(nt)==0
                            
                            u(nx,ny,nt,:)=obj.displacementsHR_stat(kx(nx),ky(ny),Pt(nx,ny,nt),ITM_PROPERTIES);
                        else
                            
                            u(nx,ny,nt,:)=obj.displacementsHR_last(kx(nx),ky(ny),omega(nt),Pt(nx,ny,nt),ITM_PROPERTIES);
                        end
                    end
                end
            end
            
            uy=(fftshift(ifft(ifftshift(u),[],2)));
            
            ux=(fftshift(ifft(ifftshift(uy),[],1)));
        end
        
        
        function pt = computeLoadSoilSuperElement(obj,Nx,Nt,dx,frequency)
            % compute load to estimate stiffness matrix
            p=zeros(Nx,Nx,Nt);
            p0=zeros(Nx,Nx);
            
            p0(Nx/2+1,Nx/2+1)=-1/(dx*dx);
            
            dt=1/(Nt*frequency);
            
            % Zeitspektrum [s]
            t=-Nt/2*dt:dt:(Nt/2-1)*dt;
            
            for st=1:Nt
                p(:,:,st)=p0*cos(2*pi*frequency*t(st));
            end
            
            % point load in wavenumber-frequency domain (kx,ky,omega);
            pt=fftshift(fft(fftshift(p),[],3));
            pt=fftshift(fft(fftshift(pt),[],2));
            pt=fftshift(fft(fftshift(pt),[],1));
            
        end
        
        
        function u=displacementsHR_stat(obj,kx,ky,Pt,ITM_PROPERTIES)
            % compute displacement for static case (omega=0)(Stiffness for each combination of kx and ky in an infinite half space
            
            mue1 = ITM_PROPERTIES(1);
            lambda1 = ITM_PROPERTIES(2);
            h1 = ITM_PROPERTIES(5);
            
            mue2 = ITM_PROPERTIES(6);
            lambda2 = ITM_PROPERTIES(7);
            
            % Substitution
            kr_2=kx^2+ky^2;
            kr=abs(sqrt(kx^2+ky^2));
            sum1=lambda1+mue1;
            sum2=lambda2+mue2;
            
            z=0;
            
            % K-Matrices
            K1o=[ 1i*kx*z                                  1i*kx*z                                  0                0             -kr              kr                   % u_x
                1i*ky*z                                  1i*ky*z                                  kr              -kr             0               0                    % u_y
                kr*z-(lambda1+3*mue1)/sum1             -kr*z-(lambda1+3*mue1)/sum1             -1i*ky            -1i*ky           1i*kx            1i*kx                 % u_z
                -(2*kx^2*z+2*lambda1*kr/sum1)           -(2*kx^2*z-2*lambda1*kr/sum1)            0                0             -2*1i*kx*kr       2*1i*kx*kr            % sigma_x
                -(2*ky^2+2*lambda1*kr/sum1)             -(2*ky^2-2*lambda1*kr/sum1)              2*1i*ky*kr       -2*1i*ky*kr      0               0                    % sigma_y
                2*kr_2*z-2*(lambda1+2*mue1)*kr/sum1     2*kr_2*z+2*(lambda1+2*mue1)*kr/sum1    -2*1i*ky*kr        2*1i*ky*kr      2*1i*kx*kr      -2*1i*kx*kr            % sigma_z
                -2*kx*ky*z                              -2*kx*ky*z                               1i*kx*kr         -1i*kx*kr       -1i*ky*kr         1i*ky*kr              % tau_xy
                2*1i*ky*(kr*z-mue1/sum1)                -2*1i*ky*(kr*z+mue1/sum1)                 kr_2+ky^2        kr_2+ky^2     -kx*ky          -kx*ky                % tau_yz
                2*1i*kx*(kr*z-mue1/sum1)                -2*1i*kx*(kr*z+mue1/sum1)                 kx*ky            kx*ky         -(kr_2+kx^2)    -(kr_2+kx^2)];        % tau_zx
            K1o(4:9,:)=K1o(4:9,:)*mue1;
            
            z=h1;
            
            K1u=[ 1i*kx*z                                  1i*kx*z                                  0                0             -kr              kr                   % u_x
                1i*ky*z                                  1i*ky*z                                  kr              -kr             0               0                    % u_y
                kr*z-(lambda1+3*mue1)/sum1             -kr*z-(lambda1+3*mue1)/sum1             -1i*ky            -1i*ky           1i*kx            1i*kx                 % u_z
                -(2*kx^2*z+2*lambda1*kr/sum1)           -(2*kx^2*z-2*lambda1*kr/sum1)            0                0             -2*1i*kx*kr       2*1i*kx*kr            % sigma_x
                -(2*ky^2+2*lambda1*kr/sum1)             -(2*ky^2-2*lambda1*kr/sum1)              2*1i*ky*kr       -2*1i*ky*kr      0               0                    % sigma_y
                2*kr_2*z-2*(lambda1+2*mue1)*kr/sum1     2*kr_2*z+2*(lambda1+2*mue1)*kr/sum1    -2*1i*ky*kr        2*1i*ky*kr      2*1i*kx*kr      -2*1i*kx*kr            % sigma_z
                -2*kx*ky*z                              -2*kx*ky*z                               1i*kx*kr         -1i*kx*kr       -1i*ky*kr         1i*ky*kr              % tau_xy
                2*1i*ky*(kr*z-mue1/sum1)                -2*1i*ky*(kr*z+mue1/sum1)                 kr_2+ky^2        kr_2+ky^2     -kx*ky          -kx*ky                % tau_yz
                2*1i*kx*(kr*z-mue1/sum1)                -2*1i*kx*(kr*z+mue1/sum1)                 kx*ky            kx*ky         -(kr_2+kx^2)    -(kr_2+kx^2)];        % tau_zx
            K1u(4:9,:)=K1u(4:9,:)*mue1;
            
            z=0;
            
            K2o=[ 1i*kx*z                                  1i*kx*z                                  0                0             -kr              kr                   % u_x
                1i*ky*z                                  1i*ky*z                                  kr              -kr             0               0                    % u_y
                kr*z-(lambda2+3*mue2)/sum2             -kr*z-(lambda2+3*mue2)/sum2             -1i*ky            -1i*ky           1i*kx            1i*kx                 % u_z
                -(2*kx^2*z+2*lambda2*kr/sum2)           -(2*kx^2*z-2*lambda2*kr/sum2)            0                0             -2*1i*kx*kr       2*1i*kx*kr            % sigma_x
                -(2*ky^2+2*lambda2*kr/sum2)             -(2*ky^2-2*lambda2*kr/sum2)              2*1i*ky*kr       -2*1i*ky*kr      0               0                    % sigma_y
                2*kr_2*z-2*(lambda2+2*mue2)*kr/sum2     2*kr_2*z+2*(lambda2+2*mue2)*kr/sum2    -2*1i*ky*kr        2*1i*ky*kr      2*1i*kx*kr      -2*1i*kx*kr            % sigma_z
                -2*kx*ky*z                              -2*kx*ky*z                               1i*kx*kr         -1i*kx*kr       -1i*ky*kr         1i*ky*kr              % tau_xy
                2*1i*ky*(kr*z-mue2/sum2)                -2*1i*ky*(kr*z+mue2/sum2)                 kr_2+ky^2        kr_2+ky^2     -kx*ky          -kx*ky                % tau_yz
                2*1i*kx*(kr*z-mue2/sum2)                -2*1i*kx*(kr*z+mue2/sum2)                 kx*ky            kx*ky         -(kr_2+kx^2)    -(kr_2+kx^2)];        % tau_zx
            K2o(4:9,:)=K2o(4:9,:)*mue2;
            
            % K-Matrix layer 1 top
            H1o=ones(9,1)*[exp(-kr*h1)  1  exp(-kr*h1)  1  exp(-kr*h1)  1];
            K1o=K1o.*H1o;
            
            % K-Matrix layer 1 bottom
            H1u=ones(9,1)*[1  exp(-kr*h1)  1  exp(-kr*h1)  1  exp(-kr*h1)];
            K1u=K1u.*H1u;
            
            % K-Matrix layer 2 top
            H2o=ones(9,1)*[1 1 1 1 1 1];
            K2o=K2o.*H2o;
            
            % system of equations
            K=[K1o([6 8 9],:)         zeros(3,3)
                K1u([1 2 3 6 8 9],:)  -K2o([1 2 3 6 8 9],[2,4,6])];
            
            % load
            Px=[0;0;Pt;0;0;0;0;0;0];
            Py=[0;Pt;0;0;0;0;0;0;0];
            Pz=[Pt;0;0;0;0;0;0;0;0];
            
            % unknowns
            Cx=linsolve(K,Px);
            Cy=linsolve(K,Py);
            Cz=linsolve(K,Pz);
            
            % FOURIER-transformed displacement (ky,z,kx,omega)
            u=zeros(9,1);
            u(1:3)=K1o(1:3,:)*Cx(1:6,:);
            u(4:6)=K1o(1:3,:)*Cy(1:6,:);
            u(7:9)=K1o(1:3,:)*Cz(1:6,:);
        end
        
        function u=displacementsHR_last(obj,kx,ky,omega,Pt,ITM_PROPERTIES)
            % compute displacement for each combination of kx,ky and omega in an infinite half space
            
            % soil parameter
            mue1 = ITM_PROPERTIES(1); lambda1 = ITM_PROPERTIES(2); cp1 = ITM_PROPERTIES(3);
            cs1 = ITM_PROPERTIES(4); h1 = ITM_PROPERTIES(5);
            
            mue2 = ITM_PROPERTIES(6);
            lambda2 = ITM_PROPERTIES(7);
            cp2 = ITM_PROPERTIES(8);
            cs2 = ITM_PROPERTIES(9);
            
            % wavenumber for P- and S-wave
            kp1=omega/cp1;
            ks1=omega/cs1;
            kp2=omega/cp2;
            ks2=omega/cs2;
            
            % Substitution
            kr_2=kx^2+ky^2;
            
            % Exponentialcoefficients
            lambda1_1=sqrt(kr_2-kp1^2);
            lambda2_1=sqrt(kr_2-ks1^2);
            lambda1_2=sqrt(kr_2-kp2^2);
            lambda2_2=sqrt(kr_2-ks2^2);
            
            % K-Matrices
            K1=[ 1i*kx                          1i*kx                          0                   0                  -lambda2_1            lambda2_1             % u_x
                1i*ky                          1i*ky                          lambda2_1           -lambda2_1            0                   0                    % u_y
                lambda1_1                     -lambda1_1                     -1i*ky               -1i*ky                1i*kx                1i*kx                 % u_z
                -(2*kx^2+lambda1/mue1*kp1^2)   -(2*kx^2+lambda1/mue1*kp1^2)    0                   0                  -2*1i*kx*lambda2_1     2*1i*kx*lambda2_1      % sigma_x
                -(2*ky^2+lambda1/mue1*kp1^2)   -(2*ky^2+lambda1/mue1*kp1^2)    2*1i*ky*lambda2_1    -2*1i*ky*lambda2_1     0                   0                    % sigma_y
                2*kr_2-ks1^2                  2*kr_2-ks1^2                 -2*1i*ky*lambda2_1     2*1i*ky*lambda2_1     2*1i*kx*lambda2_1    -2*1i*kx*lambda2_1      % sigma_z
                -2*kx*ky                      -2*kx*ky                       1i*kx*lambda2_1      -1i*kx*lambda2_1      -1i*ky*lambda2_1       1i*ky*lambda2_1        % tau_xy
                2*1i*ky*lambda1_1              -2*1i*ky*lambda1_1               lambda2_1^2+ky^2     lambda2_1^2+ky^2    -kx*ky              -kx*ky                % tau_yz
                2*1i*kx*lambda1_1              -2*1i*kx*lambda1_1               kx*ky               kx*ky              -(lambda2_1^2+kx^2)  -(lambda2_1^2+kx^2)];  % tau_zx
            K1(4:9,:)=K1(4:9,:)*mue1;
            
            K2=[ 1i*kx                          1i*kx                          0                   0                  -lambda2_2            lambda2_2             % u_x
                1i*ky                          1i*ky                          lambda2_2           -lambda2_2            0                   0                    % u_y
                lambda1_2                     -lambda1_2                     -1i*ky               -1i*ky                1i*kx                1i*kx                 % u_z
                -(2*kx^2+lambda2/mue2*kp2^2)   -(2*kx^2+lambda2/mue2*kp2^2)    0                   0                  -2*1i*kx*lambda2_2     2*1i*kx*lambda2_2      % sigma_x
                -(2*ky^2+lambda2/mue2*kp2^2)   -(2*ky^2+lambda2/mue2*kp2^2)    2*1i*ky*lambda2_2    -2*1i*ky*lambda2_2     0                   0                    % sigma_y
                2*kr_2-ks2^2                  2*kr_2-ks2^2                 -2*1i*ky*lambda2_2     2*1i*ky*lambda2_2     2*1i*kx*lambda2_2    -2*1i*kx*lambda2_2      % sigma_z
                -2*kx*ky                      -2*kx*ky                       1i*kx*lambda2_2      -1i*kx*lambda2_2      -1i*ky*lambda2_2       1i*ky*lambda2_2        % tau_xy
                2*1i*ky*lambda1_2              -2*1i*ky*lambda1_2               lambda2_2^2+ky^2     lambda2_2^2+ky^2    -kx*ky              -kx*ky                % tau_yz
                2*1i*kx*lambda1_2              -2*1i*kx*lambda1_2               kx*ky               kx*ky              -(lambda2_2^2+kx^2)  -(lambda2_2^2+kx^2)];  % tau_zx
            K2(4:9,:)=K2(4:9,:)*mue2;
            
            % K-Matrix layer 1 top
            H1o=ones(9,1)*[exp(-lambda1_1*h1)  1  exp(-lambda2_1*h1)  1  exp(-lambda2_1*h1)  1];
            K1o=K1.*H1o;
            
            % K-Matrix layer 1 bottom
            H1u=ones(9,1)*[1  exp(-lambda1_1*h1)  1  exp(-lambda2_1*h1)  1  exp(-lambda2_1*h1)];
            K1u=K1.*H1u;
            
            % K-Matrix layer 2 top
            H2o=ones(9,1)*[1 1 1 1 1 1];
            K2o=K2.*H2o;
            
            % system of equations
            K=[K1o([6 8 9],:)         zeros(3,3)
                K1u([1 2 3 6 8 9],:)  -K2o([1 2 3 6 8 9],[2,4,6])];
            
            % load
            Px=[0;0;Pt;0;0;0;0;0;0];
            Py=[0;Pt;0;0;0;0;0;0;0];
            Pz=[Pt;0;0;0;0;0;0;0;0];
            
            % unknowns
            Cx=linsolve(K,Px);
            Cy=linsolve(K,Py);
            Cz=linsolve(K,Pz);
            
            % FOURIER-transformed displacement (ky,z,kx,omega)
            u=zeros(9,1);
            u(1:3)=K1o(1:3,:)*Cx(1:6,:);
            u(4:6)=K1o(1:3,:)*Cy(1:6,:);
            u(7:9)=K1o(1:3,:)*Cz(1:6,:);
            
        end
        
    end
    
end

