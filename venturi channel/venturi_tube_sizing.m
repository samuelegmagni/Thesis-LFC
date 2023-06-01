%% Venturi tube's size 0.75 inches
L=47*1e-3;
d1 = 16.05*1e-3;
d2 = 7*1e-3;
d_conv = d1;
d_div = d1;
d_t = d2;
alpha_conv = 7;
alpha_div = 5;
l_conv1 = cotd(alpha_conv)*d_conv*0.5;
l_conv2 = cotd(alpha_conv)*d_t*0.5;
l_conv = l_conv1-l_conv2;
l_div1 = cotd(alpha_div)*d_conv*0.5;
l_div2 = cotd(alpha_div)*d_t*0.5;
l_div = l_div1 - l_div2;


%% Venturi tube's size 12 mm
L=49*1e-3;
d1 = 9*1e-3;
d2 = 5.5*1e-3;
d_conv = d1;
d_div = d1;
d_t = d2;
alpha_conv = 7;
alpha_div = 5;
l_conv1 = cotd(alpha_conv)*d_conv*0.5;
l_conv2 = cotd(alpha_conv)*d_t*0.5;
l_conv = l_conv1-l_conv2;
l_div1 = cotd(alpha_div)*d_conv*0.5;
l_div2 = cotd(alpha_div)*d_t*0.5;
l_div = l_div1 - l_div2;