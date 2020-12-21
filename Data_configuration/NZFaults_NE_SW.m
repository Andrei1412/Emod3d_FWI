function [NZfault_xy,NZfault_z]=NZFaults_NE_SW

if nargin==0
    PtSpacing=1;
end
% Gives the NZ coastline in Lat Long

data=[nan nan
28    1
29    2
30    3
30    4
31    5    
nan nan
68    1
68    2
68    3
68    4
68    5
% nan nan
% 0    1
% 1    2
% 3    3
% 5    4
% 7    5
nan nan
48    1
49    2
50    3
51    4
52    5
];

NZfault_xy=data(:,1)+1; NZfault_z=data(:,2);
%end of file
end