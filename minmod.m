function AQ=minmod(M1,M2,M3)

    if ((M1/abs(M1))==1) && ((M2/abs(M2))==1) && ((M3/abs(M3))==1)
        M=[M1 M2 M3];
        AQ=min(M);
    elseif ((M1/abs(M1))==-1) && ((M2/abs(M2))==-1) && ((M3/abs(M3))==-1)
        M=[abs(M1) abs(M2) abs(M3)];
        AQ=min(M);
    else 
        AQ=0;
    end
     
end

