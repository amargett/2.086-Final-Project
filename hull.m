%
% |The hull parameterization is defined in five chunks:|
% 
% 
% 
% |1) General Hull Form.| 
% 
% 
% 
% |This includes metrics such as the beam  at the deck, the| 
% 
% |hull taper, and the stern taper.|
% 
% 
% 
% |2) The parallel mid body cross section.| 
% 
% 
% 
% |This includes parameters for the flare, chine radius, deadrise, and keel 
% radius|
% 
% 
% 
% |3) The Bow Form.| 
% 
% 
% 
% |This includes functions that define the bow rise (rake), the profile drift 
% angle with respect to depth,|
% 
% |the profile of rocker at the front of the ship, and the profile of the location 
% where the full breadth| 
% 
% |is achieved for a given depth.|
% 
% 
% 
% |4) The Stern Form|
% 
% 
% 
% |This includes ADD DETAIL HERE WHEN FIGURED OUT|
% 
% 
% 
% |5) Bulb Forms|
% 
% 
% 
% |This includes HOW TO FORM A BULBOUS BOW AND A BULBOUS STERN|
% 
% |"""|

classdef hull
    properties
        LOA; Lb; Ls; Bd; Dd; Bs; WL; Bc; Dc; Beta; Rc; Rk; BOW; BK; Kappa; DELTA_BOW; DRIFT; Lm; KEEL; PROF;
    end
    methods
        function obj = hull(params)
            obj.LOA = params(1); 
            obj.Lb = params(2); 
            obj.Ls = params(3); 
            obj.Bd = params(4);
            obj.Dd = params(5); 
            obj.Bs = params(6); 
            obj.WL = params(7);  
            obj.Bc = params(8);  
            obj.Dc = params(9);  
            obj.Beta = params(10); 
            obj.Rc = params(11); 
            obj.Rk = params(12);
            obj.BOW = [params(13) params(14) 0];  
            obj.BK = [0 params(15) 0]; 
            obj.Kappa = params(16); 
            obj.DELTA_BOW = [params(17) params(18) 0]; 
            obj.DRIFT = [params(19) params(20) params(21)]; 
        end


%
% '''|
% 
% |=======================================================================|
% 
% |Section 2: Cross Section|
% 
% |=======================================================================|
% 
% 
% 
% |The Cross Section is defined by the following inputs:|
% 
% |0) Bd   -> The Beam at the Deck in [m] or fraction of LOA|
% 
% |1) Dd   -> The Depth of the Deck in [m] or fraction of LOA|
% 
% |2) Bc   -> The Beam at the Chine  (intersection) in [m] or fraction of LOA|
% 
% |3) Dc   -> The Depth of the Chine (intersection) in [m] or fraction of LOA|
% 
% |4) Beta -> The deadrise angle in degrees| 
% 
% |5) Rc   -> The Chine Radius in [m] or fraction of LOA|
% 
% |6) Rk   -> The keel Radius in [m] or fraction of LOA|
% 
% 
% 
% |Constraints/ NOTES to ensure realistic sizing/ shape of a hull:| 
% 
% |0) 0 <= Dc < Dd|
% 
% |1) 0 <= Beta <= 90|
% 
% |2) Rc and Rk are agebraically limited to ensure that the radius can exist 
% with the| 
% 
% |given Bd,Dd,BcdC, and Beta values.| 
% 
% 
% 
% |'''|

function halfbeam_val = halfBeam_Midbody(obj, z)
    % finds y position given z position of the parallel midbody
    % upper gunwhale line: Az + By + C = 0, where UG = [A, B, C] 
    A = [obj.Dc obj.Bc 1.0 ;
        obj.Dd obj.Bd 1.0 ;
        1.0 1.0 1.0;];
    b = [0.0; 1.0; 1.0;];
    UG = A\b;
    
    % lower gunwhale line: Az + By + C = 0 where LG = [A,B,C]. -B/A = tan(Beta)
    A = [obj.Dc obj.Bc 1.0 ; 
        sind(obj.Beta*pi/180) cosd(obj.Beta*pi/180) 0.0 ; 
        1.0 1.0 1.0; ]; 
    b = [0.0; 0.0; 1.0;]; 
    LG = A\b; 
    
    A1 = UG(1); 
    B1 = UG(2); 
    theta = atan(-B1/A1); 
    A2 = LG(1); 
    B2 = LG(2); 
    
    % find position of the center of the chine radius and intersections with
    % the chine radius
    
    Rc_Center = zeros(2); % (y,z) pair for center
    Rc_UG_int = zeros(2); % (y,z) pair for intersection of chine radius and UG line
    Rc_LG_int = zeros(2); % (y,z) pair for intersection of chine radius and LG line
    
    
    A = [0.0 0.0 B1 A1 0.0 0.0 ; 
        0.0 0.0 0.0 0.0 B2 A2 ; 
        0.0 1.0 0.0 0.0 0.0 -1.0; 
        -1.0 0.0 0.0 0.0 1.0 0.0 ; 
        0.0 1.0 0.0 -1.0 0.0 0.0 ; 
        -1.0 0.0 1.0 0.0 0.0 0.0 ; 
        ]; 
    b = [-UG(3); -LG(3); obj.Rc*cosd(obj.Beta*pi/180); 
        obj.Rc*sind(obj.Beta*pi/180); obj.Rc*cosd(theta); obj.Rc*sind(theta);]; 
    
    c = A\b; 
    
    Rc_Center = c(1:2); % center of chine radius
    Rc_UG_int  =c(3:4); % intersection btwn chine radius and upper gunwhale
    Rc_LG_int = c(5:6); % intersection of chine radius and lower gunwhale
            
    % find pos of center of the chine radius and intersections with the keel
    % radius
    
    Rk_Center = zeros(2); 
    Rk_LG_int = zeros(2); 
    Z_bottom = 0.0; 
    A = [-1.0, 0.0, 1.0, 0.0, 0.0;
          0.0, 1.0, 0.0, -1.0, 0.0;
          0.0, 0.0, B2, A2, 0.0; 
          0.0, 1.0, 0.0, 0.0, -1.0 ;
          1.0, 0.0, 0.0, 0.0, 0.0 ;]; 

    b = [obj.Rk*sind(obj.Beta*pi/180); obj.Rk*cosd(obj.Beta*pi/180); 
        -LG(3); obj.Rk*(1/2 + 1/2*sind(obj.Rk)); -obj.Rk * (1/2 - 1/2*sind(obj.Rk));]; 
    
    c = A\b; 
    
    Rk_Center = c(1:2); 
    Rk_LG_int = c(3:4); 
    Z_bottom = c(5);
    
    %%% translate Z
    Z = Z_bottom; 
    obj.Dc = obj.Dc - Z; 
    Rc_Center(2) = Rc_Center(2) - Z; 
    Rk_Center(2) = Rk_Center(2) - Z; 
    Rc_UG_int(2) = Rc_UG_int(2) - Z;
    Rc_LG_int(2) = Rc_LG_int(2) - Z; 
    Rk_LG_int(2) = Rk_LG_int(2) - Z; 
    UG(3) = -(UG(1)*obj.Dc + UG(2)*obj.Bc); 
    LG(3) = -(LG(1)*obj.Dc + LG(2)*obj.Bc); 
    Z_bottom = Z_bottom - Z; 
    
    % calculate the half beam of the cross section at a given height, z
    % If 0 > z or Dd < z, returns -1 as an error
    if (z< 0.0)
        halfbeam_val = -1;
    end
    if z> obj.Dd
        halfbeam_val =  -1; 
    end
    if (z >= 0.0) && (z<Rk_LG_int(2))
        halfbeam_val = sign(obj.Rk)*sqrt(obj.Rk^2 - (z-Rk_Center(2))^2) + Rk_Center(1); 
    end
    if (z>= Rk_LG_int(2)) && (z< Rc_LG_int(2))
        halfbeam_val = -(LG(1) * z + LG(3))/LG(2); 
    end
    if (z>= Rc_LG_int(2)) && (z < Rc_UG_int(2))
        halfbeam_val = sqrt(obj.Rc^2 - (z-Rc_Center(2))^2) + Rc_Center(1);  
    end
    if (z>= Rc_UG_int(2)) && (z<= obj.Dd)
        halfbeam_val = -(UG(1)*z + UG(3))/UG(2); 
    end
end





%
% '''|
% 
% |=======================================================================|
% 
% |Section 3: Bow Form|
% 
% |=======================================================================|
% 
% 
% 
% |The Bow Form is defined by the following inputs:|
% 
% |0) Dd   -> The Depth of the Deck in [m] or fraction of LOA|
% 
% |1) Abow  -> The z^2 term for Bow(z) that defines the profile of the bowrise|
% 
% |2) Bbow  -> The z term for Bow(z) that defines the profile of the bowrise|
% 
% |3) BK_z -> The Z Point of the intersection of the Bow rise and keel rise 
% as percentage of Dd|
% 
% |4) Kappa-> The X position where the Keel rise begins. percentage of Lb|
% 
% |5) Adel -> z^2 term for delta(z), the x position where the max Beam is achieved 
% for a given height,  x is normalized to be between 0 and 1 as a percentage of 
% Lb|
% 
% |6) Bdel -> z term for delta(z), the x position where the max Beam is achieved 
% for a given height|
% 
% |7) Adrft-> z^2 term for drift(z), the drift angle along the bowrise and keel 
% rise|
% 
% |8) Bdrft-> z term for drift(z), the drift angle along the bowrise and keel 
% rise|
% 
% |9) Cdrft-> const term for drift(z), the drift angle along the bowrise and 
% keel rise|
% 
% 
% 
% |These Parameters solve for 4 functions:|
% 
% |0) Bow(z)   -> gives the X position of the bow rise in the form Az^2 + Bz 
% + C|
% 
% |1) Keel(x)  -> gives the z height of the keel rise with respect to X in the 
% form A*(X-Kappa*Lb)^2|
% 
% |2) Delta(z) -> gives the x position between 0 and Lb where the full breadth 
% is achieved for a given z: A(z/Dd)^2 + B(z/Dd) + C = x/Lb|
% 
% |3) Drift(z) -> gives the drift angle of the bow for a given z: Az^2 + Bz 
% + C|
% 
% 
% 
% |These four functions define the following curve for each z:|
% 
% |halfBeam_Bow(x) = Y(x) = A*sin(w*x + a) + B for all z between 0 and Dd|
% 
% |Since we know two points and the derivatives of  those two points|
% 
% 
% 
% |Constraints/ NOTES to ensure realistic sizing/ shape of a hull:| 
% 
% |0) Kappa*Lb < delta(z=0)|
% 
% |1) 0 < drift(z) < 90 for 0 <= z <= Dd (only need to check at z = 0, Dd, and 
% -B/(2*A) if within range of z )|
% 
% |2) 0 <= BK_x < Kappa*Lb|
% 
% |3) 0 <= BK_z < Dd|
% 
% |4) delta(z) > Bow(z) and Keel(z) for 0 <= z <= Dd| 
% 
% 
% 
% |'''|

    function pos = bowrise(obj,z)
        % returns the x position of the bowrise for a given z for BK_z <= z <= Dd
        pos = obj.BOW(1)*z^2 + obj.BOW(2)*z + obj.BOW(2); 
    end
    
    function pos = keelrise_bow(obj,z)
    %returns the x position of the keelrise  at the bow for a given z for 0 <= z <= Bk_z
        pos = -sqrt(z/obj.KEEL) + obj.Kappa*obj.Lb; 
    end
    
    function pos = delta_bow(obj,z)
        %returns the x position where the full cross section width is achieved for a given z for 0 <= z <= Dd
        pos = obj.Lb + obj.DELTA_BOW(1)*z^2 + obj.DELTA_BOW(2)*z + obj.DELTA_BOW(3); 
    end
    
    function val = drift(obj,z)
        val = obj.DRIFT(1) * z^2 + obj.DRIFT(2)*z + obj.DRIFT(3); 
    end
    
    function profile = bow_profile(obj,z)
        if z <= obj.BK(2)
            X1 = keelrise_bow(obj,z); 
        else
            X1 = bowrise(obj, z); 
        end
        profile = X1; 
       end
    
    function soln = solve_waterline_bow(obj,z)
    %this function solves for a cubic function: y(half beam) = Ax^3 + Bx^2 + CX + D for the half beam of the profile between the bow/keel rise and delta for a given z for 0 <= z <= D
        X1 = bow_profile(obj, z); 
        X2 = delta_bow(obj, z); 
        Y2 = halfBeam_Midbody(obj, z); 
    
        A = [X1^3, X1^2, X1, 1; 
            3*X1^2, 2*X1, 1, 0; 
            X2^3, X2^2, X2, 1; 
            3*X2^2, 2*X2, 1, 0;]; 
        b = [0; tand(pi*drift(obj, z)); Y2; 0;]; 
        soln = A\b; 
    end
        
    function halfBeam = halfBeam_Bow(obj,x)
    %#returns the halfbeam along the bow taper between the bow/keel rise and delta(z), PROF is the output of solve)waterline_bow(z)
        halfBeam = obj.PROF(1)*x^3 + obj.PROF(2)*x^2 + obj.PROF(3)*x + obj.PROF(4); 
    end
    
    function y = waterline_bow(obj, x, z)
        % finds points along waterline of bow
        obj.BK(2) = obj.BK(2) * obj.Dd; 
        if obj.BOW(1) == 0 
            Zv = -1.0; 
        else 
            Zv = obj.BOW(2)/(2*obj.BOW(1));
        end
    
        C = [obj.BOW(1)*obj.Dd^2 + obj.BOW(1)*obj.Dd
            obj.BOW(1) + obj.BK(2)^2 + obj.BOW(2) 
            obj.BOW(1)*Zv^2 + obj.BOW(2)*Zv]; 

        if (Zv >= obj.BK(2)*obj.Dd) && (Zv <= obj.Dd)
            obj.BOW(3) = -min(C); 
        else
            obj.BOW(2) = -min(C());
        end
    
    
        %Calculate the C for the Delta equation, where C is the constant such that max(Delta(z)) = 0 between 0 and Dd
        if obj.DELTA_BOW(1) == 0
            Zv = -1.0; 
        else 
            Zv = -obj.DELTA_BOW(2)/(2*obj.DELTA_BOW(1)); % Find Z of vertex of Delta(z)

            C = [(obj.DELTA_BOW(1)*obj.Dd^2 + obj.DELTA_BOW(2)*obj.Dd);
            0.0; (obj.DELTA_BOW(1) * Zv^2 + obj.DELTA_BOW(2) * Zv);];
        end
        if (Zv >= 0.0) && (Zv <= obj.Dd)
            obj.DELTA_BOW(3) = -max(C); 
        else
            obj.DELTA_BOW(3) = -max(C(1:2)); 
        end
    
    
    %  This generates a y value to detail the curvature of the bow taper for a given z, for 0 <= z <= D

        obj.PROF = solve_waterline_bow(obj, z); 
        y = halfBeam_Bow(obj, x);

    end
    end
end


