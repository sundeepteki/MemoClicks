function [p,anovatab,stats] = st_anova1(x,y,z)

%{
   Sundeep Teki. 08.07.15
   function to perform 1 way anova for 3 inputs
   usage: [a,b,c] = st_anova1(t1,t2,t3)

%}


t = [x;y;z]';
[p,anovatab,stats] = anova1(t,[],'off');

end
