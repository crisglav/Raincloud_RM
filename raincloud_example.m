% Example of raincloud plots
f = figure;
n = 250;
d{1} = (exprnd(5, 1, n) + 15)';
d{2} = ((randn(1, n) *5) + 20)';
d{3} = ((randn(1, n) *5) + 20)';

cl = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4];
h   = rm_raincloud_cg(d', 'colours',cl,'box_on',1, 'box_dodge',1,'line_width',1);

