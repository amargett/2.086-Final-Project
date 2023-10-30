classdef boat1
    properties
        LOA; Lb; Ls; Bd; Dd; Bc; Dc; points;
    end
    methods
        function obj = boat1()
        obj.LOA = 10.3;
        obj.Lb = 3.6;
        obj.Ls = 1.8;
        obj.Bd = 2.73;
        obj.Dd = 2.14;
        obj.Bc = 2.37;
        obj.Dc = 9.67;
        obj.points = [0.9, -1.25, 1.65; 6.82, -0.02, 1.19; 4.84, 0.01, 1.37; 4.96, 0.03, 1.81; 6.36, 0.0, 1.52; 0.58, -1.35, 1.87; 7.57, -0.01, 1.63; 5.01, 0.01, 1.4; 3.85, 0.01, 1.89; 7.0, 0.01, 2.03; 9.89, -0.07, 1.39; 1.04, -1.1, 1.64; 8.73, -0.06, 1.26; 9.08, -0.06, 1.4; 9.37, -0.07, 1.23; 7.46, -0.01, 1.75; 7.4, -0.01, 1.72; 5.25, 0.0, 1.39; 5.91, 0.0, 1.52; 5.69, 0.02, 1.75; 4.87, 0.03, 1.8; 7.52, -0.03, 1.35; 6.76, 0.0, 1.59; 4.11, -0.01, 1.48; 3.77, -0.05, 1.36; 2.13, -0.43, 1.53; 3.97, 0.04, 2.13; 8.43, -0.01, 1.89; 8.15, -0.03, 1.54; 7.34, -0.02, 1.44; 4.8, 0.01, 1.49; 6.04, -0.01, 1.22; 5.07, 0.0, 1.23; 4.85, 0.0, 1.27; 5.33, 0.0, 1.27; 7.07, -0.03, 1.15; 3.38, -0.03, 1.83; 
            7.78, -0.61, 0.99; 8.96, -1.22, 0.75; 8.29, -0.61, 0.99; 8.63, -0.79, 0.92; 8.87, -0.42, 1.06; 7.6, -1.22, 0.75; 8.38, -0.29, 1.11; 8.92, -0.28, 1.12; 9.59, -0.29, 1.11; 8.35, -0.31, 1.1; 7.54, -0.77, 0.93; 7.54, -1.19, 0.76; 9.27, -1.16, 0.77; 8.23, -0.29, 1.11; 8.7, -0.91, 0.87; 8.04, -0.79, 0.92; 7.85, -0.73, 0.94; 8.61, -1.16, 0.77; 8.77, -1.28, 0.72; 8.9, -0.78, 0.92; 7.43, -0.35, 1.09; 7.47, -1.09, 0.8; 8.84, -0.74, 0.94; 7.92, -0.72, 0.95; 9.0, -0.44, 1.05; 9.37, -1.23, 0.74; 7.36, -0.9, 0.88; 8.12, -0.3, 1.11; 7.39, -0.37, 1.08; 8.29, -0.39, 1.07; 9.32, -1.17, 0.77; 8.68, -0.74, 0.94; 7.38, -0.45, 1.05; 9.1, -0.77, 0.92; 8.69, -0.75, 0.93; 9.11, -0.86, 0.89; 7.69, -1.14, 0.78; 8.34, -0.14, 1.16;
            7.25, -1.07, 0.76; 6.66, -0.36, 1.04; 6.59, -0.05, 1.15; 6.55, -0.59, 0.96; 6.98, -0.15, 1.1; 6.83, -1.07, 0.76; 6.29, -0.2, 1.11; 7.14, -0.8, 0.86; 7.27, -1.36, 0.64; 6.66, -1.33, 0.66; 6.58, -1.11, 0.75; 6.26, -1.24, 0.71; 6.87, -1.3, 0.67; 6.32, -0.09, 1.15; 6.9, -1.19, 0.71; 6.6, -0.23, 1.09; 6.12, -0.4, 1.04; 6.99, -1.08, 0.76; 6.74, -1.1, 0.75; 6.27, -1.31, 0.68; 6.26, -0.9, 0.85; 7.25, -1.24, 0.69; 6.0, -1.1, 0.77; 6.81, -0.39, 1.02; 6.21, -1.14, 0.75; 7.28, -0.94, 0.81; 6.77, -0.84, 0.86; 6.3, -0.93, 0.83; 6.6, -1.37, 0.64; 6.37, -0.06, 1.16; 6.03, -1.33, 0.67; 6.5, -0.17, 1.11; 6.3, -0.12, 1.14; 6.89, -1.27, 0.68; 6.61, -0.44, 1.01; 6.66, -0.55, 0.97;
            2.96, -1.3, 0.75; 5.27, -0.79, 0.86; 2.41, -0.67, 1.26; 1.3, -1.0, 1.53; 4.86, -1.12, 0.72; 3.29, -1.13, 0.81; 5.04, -0.4, 1.02; 3.46, -0.4, 1.16; 5.07, -0.58, 0.95; 5.04, -0.73, 0.89; 4.03, -0.88, 0.87; 4.44, -0.31, 1.09; 3.6, -0.66, 1.02; 2.63, -0.83, 1.11; 1.36, -1.02, 1.49; 4.3, -1.28, 0.66; 4.39, -0.75, 0.9; 4.96, -1.28, 0.65; 1.95, -1.29, 1.02; 1.28, -0.96, 1.57; 2.04, -0.71, 1.36; 3.54, -1.24, 0.72; 2.17, -0.58, 1.39; 1.4, -1.01, 1.47; 5.19, -0.7, 0.9; 3.87, -1.13, 0.75; 4.43, -0.86, 0.86; 5.61, -1.21, 0.67; 5.01, -0.75, 0.88; 3.22, -0.22, 1.29; 3.53, -0.3, 1.2; 4.68, -0.09, 1.16; 1.03, -1.25, 1.55; ];
        
       obj.points(:, 1) = obj.points(:,1) - min(obj.points(:,1)); 
       obj.points(:,2) = obj.points(:,2) - min(obj.points(:,2));
       obj.points(:,3) = obj.points(:,3) - min(obj.points(:,3)); 

        end
    end
end