classdef FM_class
   %% ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    properties
        d_depth; % vertical resolution [m]
        d_alpha_mdl; % horizontal resolution [°]
        d_alpha_synt; % horizontal resolution [°]
        t_layer; % depth of the horizontal layer boundary
        theta; % fabric orientation (abgle between v1 and TR)
        AzOfst; % Offset azimuth for theta
        lambda_1; % smallest horizontal eigenvalue
        lambda_2; % largest horizontal eigenvalue
        HA; % Horizontal anisotropy (lambda2 - lambda1)
        r; % reflection ratio as a ratio (1 = no anisotropic scattering)
        rdB; % reflection ratio in dB (0 = no anisotropic scattering)
        f; % radar center frequency [Hz]
        epsX; % Dielectric permitivity in X
        epsY; % Dielectric permitivity in Y
        GammaX; % Complex amplitudes reflection coefficients X
        GammaY; % Complex amplitudes reflection coefficients Y
        delX; % Conductivity X
        delY; % Conductivity Y
        DelEpsP; % Delta Epsilon Prime -- > Dielectric anisotropy  
        c; % Speed of light [m/s]
        epsilon0; % Dielectric permitivity in a vacuum [F/m]
        P1; % First coefficient matrix
        P2; % Second coefficient matrix
        Zmx; % Maximum depth
        Z; % Depth vector
        nL; % number of layers
        omega; % Angular frequency [rad/s]
        lambda0; % Wavelength in a vacuum [m]
        k0; % [rad/m]
        mu0; % Magnetic permeability in a vacuum [Henry/m]
        j = sqrt(-1); % Complex number
        SIGNAL; % Modeled Signal
        HH;
        VH;
        HV;
        VV;
        Dta;
    end
   %% ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    methods(Static)
       %% ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       function obj = CalculateSecondaryParameters(obj)
            obj.nL = size(obj.t_layer,1); % number of layers
            %---------------------------------------------------------
            obj.Zmx = max(obj.t_layer); % Maximum depth
            obj.Z = (0:obj.d_depth:obj.Zmx)'; % depth vector
            obj.Z(1) = 1e-20;
            %---------------------------------------------------------
            obj.r = 10 .^ (obj.rdB / 20); % dB to ratio
            obj.epsX = repmat(3.15,obj.nL,1);
            obj.epsY = (obj.HA .* obj.DelEpsP) + obj.epsX;
            obj.GammaX = repmat(1e-12,obj.nL,1); 
            obj.GammaY = obj.GammaX .* obj.r;
            %---------------------------------------------------------
            obj.omega = 2 * pi * obj.f; % angular frequency [rad/s]
            obj.lambda0 = obj.c / obj.f; % wavelength in a vacuum [m]
            obj.k0 = 2 * pi / obj.lambda0; % [rad/m]
            %---------------------------------------------------------
            obj.P1 = repmat([1 0 ; 0 1],1,1,length(obj.d_alpha_mdl)); % First coefficient matrix
            obj.P2 = repmat([1 0 ; 0 1],1,1,length(obj.d_alpha_mdl)); % Second coefficient matrix 
       end
       %% ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       function obj = RunTheModel(obj)
            %---------------------------------------------------------
            iZ0 = 1; % index of the depth at the surface
            for i = 1:obj.nL
                if obj.nL == 1
                    iZ1 = length(obj.Z);
                else
                    [~,iZ1] = min(abs(obj.t_layer(i)-obj.Z)); % index of current layer boundary
                end
                %---------------------------------------------------------
                Z = obj.Z(iZ0:iZ1); % Depth vector of the current layer
                epsX = obj.epsX(i);
                epsY = obj.epsY(i);
                GammaX = obj.GammaX(i);
                GammaY = obj.GammaY(i);
                theta = obj.theta(i);
                %---------------------------------------------------------
                rdpX = epsX - (obj.j .* (obj.delX ./(obj.epsilon0 .* obj.omega)));
                kx = sqrt((obj.epsilon0 .* obj.mu0 .* rdpX .* (obj.omega.^2)) + (obj.j .* obj.mu0 .* obj.delX .* obj.omega));
                Tx = exp( (-obj.j .* obj.k0 .* obj.d_depth) + (obj.j .* kx .* obj.d_depth) );
                %---------------------------------------------------------
                rdpY = epsY - (obj.j .* (obj.delY ./(obj.epsilon0 .* obj.omega)));
                ky = sqrt((obj.epsilon0 .* obj.mu0 .* rdpY .* (obj.omega.^2)) + (obj.j .* obj.mu0 .* obj.delY .* obj.omega));
                Ty = exp( (-obj.j .* obj.k0 .* obj.d_depth) + (obj.j .* ky .* obj.d_depth) );
                %---------------------------------------------------------
                T = [Tx complex(0) ; complex(0) Ty];
                G = [GammaX 0 ; 0 GammaY];
                ti = obj.d_alpha_mdl - theta;
                sL = length(Z);
                %---------------------------------------------------------
                HH = nan(sL,size(ti,2));
                VV = nan(sL,size(ti,2));
                HV = nan(sL,size(ti,2));
                VH = nan(sL,size(ti,2));
                Signal = nan(sL,size(ti,2),4);
                %---------------------------------------------------------
                for m = 1:size(ti,2)
                    P1 = obj.P1(:,:,m);
                    P2 = obj.P2(:,:,m);
                    dg = ti(1,m);
                    %---------------------------------------------------------
                    n=(1:sL)';
                    TT = [T(1,1) T(2,2)].^n;
                    D = (exp(obj.j.*obj.k0.*Z)./(4.*pi.*Z)).^2;
                    csd = cosd(dg);
                    snd = sind(dg);
                    PRGRP = P1 * FM_class.Rmat(dg) * G * FM_class.Rmat(-dg) * P2;
                    A1 = PRGRP(1,1);
                    A2 = PRGRP(1,2);
                    A3 = PRGRP(2,1);
                    A4 = PRGRP(2,2);
                    a1 = ((csd.^2) .* TT(:,1) ) + ((snd.^2) .* TT(:,2) );
                    a2 = csd .* snd .* (TT(:,1) - TT(:,2));
                    a4 = ((snd.^2) .* TT(:,1) ) + ((csd.^2) .* TT(:,2) );
                    %---------------------------------------------------------
                    HH(:,m) = D.* ( a1.*(a1.*A1 + a2.*A3) + a2.*(a1.*A2 + a2.*A4));
                    VH(:,m) = D.* ( a2.*(a1.*A1 + a2.*A3) + a4.*(a1.*A2 + a2.*A4));
                    HV(:,m) = D.* ( a1.*(a2.*A1 + a4.*A3) + a2.*(a2.*A2 + a4.*A4));
                    VV(:,m) = D.* ( a2.*(a2.*A1 + a4.*A3) + a4.*(a2.*A2 + a4.*A4));
                    %---------------------------------------------------------
                    obj.P1(:,:,m) =  ((FM_class.Rmat(dg) * T * FM_class.Rmat(-dg)) ^ (sL)) * P1;
                    obj.P2(:,:,m) =  P2 * ((FM_class.Rmat(dg) * T * FM_class.Rmat(-dg)) ^ (sL));
                end
                %---------------------------------------------------------
                Signal(:,:,1) = HH;
                Signal(:,:,2) = VV;
                Signal(:,:,3) = HV;
                Signal(:,:,4) = VH;
                obj.SIGNAL(iZ0:iZ1,:,:) = Signal;
                iZ0 = iZ1 + 1;
            end 
       end
       %% ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       function R = Rmat(deg)
            R = [cosd(deg) -sind(deg);
                 sind(deg)  cosd(deg)];
       end
       %% ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       function obj = AzimuthSynthesizer(obj)
            HH0 = obj.SIGNAL(:,:,1);    
            VV0 = obj.SIGNAL(:,:,2);   
            HV0 = obj.SIGNAL(:,:,3);  
            VH0 = obj.SIGNAL(:,:,4);
            if size(HH0,2) == length(obj.d_alpha_synt)
                obj.HH = HH0; 
                obj.VV = VV0; 
                obj.HV = HV0; 
                obj.VH = VH0;
            elseif size(HH0,2) == 1
                hh = HH0(:,1);
                vh = VH0(:,1);
                hv = HV0(:,1);
                vv = VV0(:,1);
                %---------------------------------------------------------
                obj.HH = zeros(size(hh,1),size(obj.d_alpha_synt,2));
                obj.VV = zeros(size(obj.HH));
                obj.HV = zeros(size(obj.HH));
                obj.VH = zeros(size(obj.HH));
                %---------------------------------------------------------
                vh = hv;
                %---------------------------------------------------------
                theta_prime = obj.d_alpha_synt+obj.AzOfst;
                for kk=1:length(theta_prime)
                    t= theta_prime(kk)*pi/180;
                    obj.HH(:,kk) = (hh*cos(t)^2) + (vv*sin(t)^2) - (sin(t)*cos(t)*(hv+vh));
                    obj.VH(:,kk) = (cos(t)^2*vh) - (sin(t)^2*hv) + (sin(t)*cos(t)*(hh-vv));
                    obj.HV(:,kk) = (cos(t)^2*hv) - (sin(t)^2*vh) + (sin(t)*cos(t)*(hh-vv));
                    obj.VV(:,kk) = (vv*cos(t)^2) + (hh*sin(t)^2) + (sin(t)*cos(t)*(hv+vh));
                end  
            end
       end
       %% ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        function obj = Signal2Param(obj) 
            HH = obj.HH;
            VH = obj.VH;
            HV = obj.HV;
            VV = obj.VV;
            Z = obj.Z;
            ao = obj.d_alpha_synt;
            f = obj.f;
            dZ = mean(diff(Z));
            DenoisingFlag=[]; C_DepthWin=dZ; C_ConvWin=dZ; dtatyp = "model";
            %--------------------------------------------------------- Power anomaly & lateral phase difference
            [PA_HH,PHI_HH] = FM_class.AmpPhs(HH,ao);
            [PA_VH,PHI_VH] = FM_class.AmpPhs(VH,ao);
            [PA_HV,PHI_HV] = FM_class.AmpPhs(HV,ao);
            [PA_VV,PHI_VV] = FM_class.AmpPhs(VV,ao);
            %--------------------------------------------------------- Polarimetric Coherence
            C_HHVV = FM_class.Chhvv(HH,VV,C_DepthWin,dZ,dtatyp);
            absC = abs(C_HHVV);
            argC = angle(C_HHVV);
            reC = real(C_HHVV);
            imC = imag(C_HHVV);
            if isempty(C_ConvWin)
                C_ConvWin = C_DepthWin;
            end
            %--------------------------------------------------------- Polarimetric Coherence depth gradient
            gradC = FM_class.Chhvv_Z_derivative(C_HHVV,C_ConvWin,dZ);
            Psi = (2*299792458*sqrt(3.15)*gradC) ./ (4*pi*f*0.034);
            %--------------------------------------------------------- Output structure
            % PA,PD : HH VH HV VV
            % Coh: abs, arg, re, im , grad, Psi
            Sig.HH = HH; Sig.VH = VH; Sig.HV = HV; Sig.VV = VV;
            PA.HH = PA_HH; PA.VH = PA_VH; PA.HV = PA_HV; PA.VV = PA_VV;
            PD.HH = PHI_HH; PD.VH = PHI_VH; PD.HV = PHI_HV; PD.VV = PHI_VV;
            C.absC = absC; C.argC = argC; C.reC = reC; C.imC = imC;
            C.gradC = gradC; C.Psi = Psi;
            %--------------------------------------------------------- Denoising
%             [PA,PD] = CLASS_Denoising.RunDenoising(DenoisingFlag,PA,PD,dZ);
%             PA_HH = PA.HH; PA_VH = PA.VH; PA_HV = PA.HV; PA_VV = PA.VV;
%             PHI_HH = PD.HH; PHI_VH = PD.VH; PHI_HV = PD.HV; PHI_VV = PD.VV;
            %--------------------------------------------------------- Store
            obj.Dta = {  HH,VH,HV,VV,... % 1 2 3 4
                     PA_HH,PA_VH,PA_HV,PA_VV,... % 5 6 7 8
                     PHI_HH,PHI_VH,PHI_HV,PHI_VV,... % 9 10 11 12
                     absC,argC,reC,imC,gradC,Psi}; % 13 14 15 16 17 18
        end 
       %% ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        % Power Anomaly & Lateral Phase Difference
        function [delta_P,dphi_dtheta] = AmpPhs(s,ao)
            delta_P = 20.*log10(abs(s)./mean(abs([s]),2,'omitnan')); % Power anomaly
            delta_P(isinf(delta_P)) = nan;
            s = s.';
            PhaseDiff = nan(size(s,1)-1,size(s,2));
            for i=1:length(ao)-1
                 PhaseDiff(i,:)=angle(s(i+1,:).*conj(s(i,:)));
            end
            xdiffphase=cos(PhaseDiff');
            ydiffphase=sin(PhaseDiff'); 
            dphi_dtheta=atan2(ydiffphase,xdiffphase);
        end
       %% ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        function Chhvv = Chhvv(Shh,Svv,DptAvgC,dZ,dtatype)
            stwnd = DptAvgC/dZ;
            Chhvv1 = movsum(Shh.*conj(Svv),[stwnd],1,'omitnan');
            Chhvv2 = sqrt(movsum(abs(Shh).^2,[stwnd],1,'omitnan'));
            Chhvv3 = sqrt(movsum(abs(Svv).^2,[stwnd],1,'omitnan'));
            Chhvv = Chhvv1 ./ (Chhvv2 .* Chhvv3);
            if dtatype == "radar"
                Chhvv = conj(Chhvv); % assumption from Jordan (converting the synthetic received seignal to de-ramped signal)
            end
        end 
       %% ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        function dphi_dz = Chhvv_Z_derivative(Chhvv,ConvWin,dZ)
            ConvWin =round(ConvWin/dZ,0);
            win=gausswin(ConvWin);
            win=win/sum(win);
            R = real(Chhvv);
            I = imag(Chhvv); 
            Rz=[];
            Iz=[];
            Rlow=[];
            Ilow=[];
            for j=1:size(Chhvv,2)
                Rz(:,j)=diff(conv(R(:,j),win,'same'))/dZ;
                Iz(:,j)=diff(conv(I(:,j),win,'same'))/dZ;
                Rlow(:,j)=conv(R(:,j),win,'same');
                Ilow(:,j)=conv(I(:,j),win,'same');
            end
            Iz(size(I,1),:)=NaN;
            Rz(size(R,1),:)=NaN;
            dphi_dz =((R.*Iz-I.*Rz)./(Rlow.^2+Ilow.^2));
        end
       %% ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       function f = SetFigureSize(ss1,ss2,w,h)
            f = figure;
            set(f,'Color',[1 1 1]);
            set(f, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]); % full screen figure
            set(f, 'Units', 'centimeters');
            scrn = get(f, 'OuterPosition'); % get the size of the screen in CM
        
            wdt = scrn(3) * w;
            hgt = scrn(4) * h;
        
            s1 = scrn(3) * ss1;
            s2 = scrn(4) * ss2;
        
            set(f, 'OuterPosition', [s1, s2, wdt, hgt]); % change the figure size to the new size
            set(f, 'Units', 'Normalized');
        end
    end
end