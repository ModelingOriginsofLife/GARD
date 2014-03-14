function J = tgscolors(type)
%function J = tgscolors(type)
%some different colormaps according to Doron's whims.
%inputs:
%type - currently may be 'bone' for the beta matrix or 'jet' for the reduced_jet carpets
%outputs:
%J - returns the colormap - so then type colormap(J)
%Or with no output arguments automatically does this.
J=[];

switch lower(type)
    case {'bone'}
        J=colormap(bone);J=flipud(J);
    case {'jet'}
        J =  [  0       0    0.8971
            0    0.0057    0.9449
            0    0.0241    0.9800
            0    0.0621    0.9955
            0    0.1117    0.9995
            0    0.1647    1.0000
            0    0.2182    1.0000
            0    0.2717    1.0000
            0    0.3252    1.0000
            0    0.3787    1.0000
            0    0.4323    1.0000
            0    0.4858    1.0000
            0    0.5393    1.0000
            0    0.5928    1.0000
            0    0.6464    1.0000
            0    0.6999    1.0000
            0    0.7534    1.0000
            0    0.8069    1.0000
            0.0000    0.8604    1.0000
            0.0008    0.9132    0.9992
            0.0107    0.9568    0.9893
            0.0367    0.9843    0.9633
            0.0780    0.9965    0.9220
            0.1284    0.9997    0.8716
            0.1816    1.0000    0.8184
            0.2351    1.0000    0.7649
            0.2886    1.0000    0.7114
            0.3421    1.0000    0.6579
            0.3957    1.0000    0.6043
            0.4492    1.0000    0.5508
            0.5027    1.0000    0.4973
            0.5562    1.0000    0.4438
            0.6097    1.0000    0.3903
            0.6633    1.0000    0.3367
            0.7168    1.0000    0.2832
            0.7703    1.0000    0.2297
            0.8238    1.0000    0.1762
            0.8770    0.9997    0.1230
            0.9276    0.9967    0.0724
            0.9683    0.9839    0.0317
            0.9912    0.9533    0.0088
            0.9986    0.9072    0.0014
            0.9999    0.8550    0.0001
            1.0000    0.8015         0
            1.0000    0.7480         0
            1.0000    0.6945         0
            1.0000    0.6410         0
            1.0000    0.5874         0
            1.0000    0.5339         0
            1.0000    0.4804         0
            1.0000    0.4269         0
            1.0000    0.3734         0
            1.0000    0.3198         0
            1.0000    0.2663         0
            1.0000    0.2128         0
            0.9999    0.1593         0];

    otherwise
        disp('Unknown method.')
end

if nargout==0; colormap(J); end