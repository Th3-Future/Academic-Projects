function [BB,pqr] = Sensor(BB,pqr)

for idx = 1:3
    
 %%%Update Rate
nextSensorUpdate = 1;
   %%%our sensor paramameters
%%%Bias and Noise
MagscaleBias = (4e-7); %%T
MagFieldBias = MagscaleBias*(2*rand()-1); %%% -1 to 1 (0 to 1) (0 to 2
MagscaleNoise = (1e-5); %%T
MagFieldNoise = MagscaleNoise*(2*rand()-1); %%% -1 to 1 (0 to 1) (0 to 2
AngscaleBias = 0.01; %%rad/s
AngFieldBias = AngscaleBias*(2*rand()-1); %%% -1 to 1 (0 to 1) (0 to 2
AngscaleNoise = 0.001; %%rad/s
AngFieldNoise = AngscaleNoise*(2*rand()-1); %%% -1 to 1 (0 to 1) (0 to 2
    %%%Pollute the data
    BB(idx) = BB(idx) + MagFieldBias + MagFieldNoise;
    pqr(idx) = pqr(idx) + AngFieldBias + AngFieldNoise;
end


