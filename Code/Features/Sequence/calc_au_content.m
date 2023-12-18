function au_content = calc_au_content(RNA,RNA_start,seed_length,seed_type)
%Calculating weighted AU content

%RNA = The whole RNA strand
%RNA_start = Start position of the seed in the RNA
%seed_length = Length between 6 and 7

load('au_params.mat','au_params');

switch seed_type
    case '6mer'
        params = au_params.six;
    case '6mer-A1'
        params = au_params.six;
    case 'offset-6mer'
        params = au_params.six;
    case '7mer-A1'
        params = au_params.sevenA1;
    case 'offset-7mer'
        params = au_params.sevenA1;
    case '7mer-m8'
        params = au_params.sevenm8;
    case '8mer'
        params = au_params.eight;
end

%Parameters:
%weights_up/down = Weights for up/downstream strand
%window - length of window to check for au content
%up/down_shift - shift from beginning of miRNA (au content is not checked
%in the seed)
window = params.window_length;
up_shift = params.up_shift;
down_shift = params.down_shift;

weights_up = params.weights_up;
weights_down = params.weights_down;

RNA = lower(RNA);
RNA(RNA == 't') = 'u';

seq_up = RNA(max(1,RNA_start+seed_length+1-up_shift-window):RNA_start+seed_length-up_shift);
seq_down = RNA(RNA_start+seed_length+down_shift:min(length(RNA),RNA_start+seed_length+window+down_shift-1));

wup = weights_up(length(weights_up)-length(seq_up)+1:end);
wdown = weights_down(1:length(seq_down));

au_up_index = sort([strfind(seq_up,'a'),strfind(seq_up,'u')]);
au_down_index = sort([strfind(seq_down,'a'),strfind(seq_down,'u')]);

au_up = sum(1./wup(au_up_index));
au_down = sum(1./wdown(au_down_index));

au_content = (au_up+au_down)/(sum(1./wup)+sum(1./wdown));

end

