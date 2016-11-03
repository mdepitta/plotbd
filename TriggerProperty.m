function propval = TriggerProperty(index,defval,altval,index_alt)
%
% Utility that allows to trigger different values depending on
% 'index'. It alternates between defval and altval, returning altval only 
% for integer numbers given in vector index_alt
%
% Maurizio De Pitta', Tel Aviv, January 30th, 2012.

if intersect(index,index_alt)==index
    propval = altval;
else
    propval = defval;
end
