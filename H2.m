function [h]=H2(x)
h=[];
h=-x.*log2(x)-(1-x).*log2(1-x);
h(x==0)=0;
end