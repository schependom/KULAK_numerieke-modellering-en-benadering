function [W2,W3,W4,b2,b3,b4] = lees_parameters(c)
     W2 = zeros(2,2);
     W3 = zeros(3,2);
     W4 = zeros(2,3);
     W2(:) = c(1:4);
     W3(:) = c(5:10);
     W4(:) = c(11:16);
     b2 = c(17:18);
     b3 = c(19:21);
     b4 = c(22:23);
end