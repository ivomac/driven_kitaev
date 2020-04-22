function A = C(varargin)

m = numel(varargin);

A = C1(varargin{end-1},varargin{end});

for n = 2:m-1
    A = C1(varargin{end-n},A);
end

end

function A = C1(M0,M1)

    A = M0*M1-M1*M0;

end