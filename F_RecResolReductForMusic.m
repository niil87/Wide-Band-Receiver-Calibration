function doa = F_RecResolReductForMusic(V,D,timeStep)

    doa = 0;
    decPlaces = 15;
    while(doa == 0 && decPlaces > 5)
        V = round(V,decPlaces);
        CovMat = V*D*V';
        try
            doa = musicdoa(CovMat,1,'ScanAngles',[-90:0.1:90]);
            disp(["DOA Attempt successful for TimeStep:" + timeStep + " doasSage value:" + doa]);
        catch
            decPlaces = decPlaces - 1;
            disp(["DOA Attempt failed for TimeStep:" + timeStep + ", decreasing resolution to " + decPlaces + " decimal places"]);
            doa = 0;
        end
    end

end