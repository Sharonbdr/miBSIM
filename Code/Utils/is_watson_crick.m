function is_wc = is_watson_crick(str1,str2,is_wobble)
%Return logical array - 1 if the nts are a Watson-Crick RNA pair, 0 otherwise
%Indices of the array the same as strand1

if nargin < 3
    is_wobble = 0;
end

%Strands have to be of same length
if length(str1) ~= length(str2)
    error('Strands have to be of same length');
end

str1 = upper(str1);
str2 = upper(str2);

str2 = seqrcomplement(str2);
str1(str1 == 'T') = 'U';
str2(str2 == 'T') = 'U';

if is_wobble == 1
    is_wc = ( (str1 == str2) | (str1 == 'G' & str2 == 'A') | (str1 == 'U' & str2 == 'C') );
else
    is_wc = str1 == str2;
end

end