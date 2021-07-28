% calculation accuracy/ error
% any deviation is considered a relevant deviation regardles whether plus
% or minus

function A= fam_error(i_total,i_ideal)
A=((i_total-i_ideal)/i_ideal)*100;
                if A <=0        % any deviation is relevant
                    A=A*(-1);
                end
end
