function [ eqs ] = equationMatrixMaker( hammy,uni,t )
%This function creates a matrix/cell array of all the coupled differential
%equations to be solved with ode45

eqs = zeros(256,1);

% for a = 1:1:256
%     eqs{a,1} = @(t)uni(a);
% end

H_matrix = hammy(t);

for a = 1:1:256
    if mod(a,16) >= 1
        row = mod(a,16);
    else
        row = 16;
    end
    col = floor(((a)/16)-0.0001)+1;
    if mod((a),16) >= 1
        moda = mod((a),16);
    else
        moda = 16;
    end
    
%     eqs{a,1} = @(t) -i*matrixElements(hammy,row,t)*[uni(1+(a)-moda);...
%                         uni(2+(a)-moda);uni(3+(a)-moda);uni(4+(a)-moda);...
%                         uni(5+(a)-moda);uni(6+(a)-moda);uni(7+(a)-moda);...
%                         uni(8+(a)-moda);uni(9+(a)-moda);uni(10+(a)-moda);...
%                         uni(11+(a)-moda);uni(12+(a)-moda);uni(13+(a)-moda);...
%                         uni(14+(a)-moda);uni(15+(a)-moda);uni(16+(a)-moda)];

    eqs(a,1) = -i*H_matrix(row,:)*[uni(1+(a)-moda);...
                    uni(2+(a)-moda);uni(3+(a)-moda);uni(4+(a)-moda);...
                    uni(5+(a)-moda);uni(6+(a)-moda);uni(7+(a)-moda);...
                    uni(8+(a)-moda);uni(9+(a)-moda);uni(10+(a)-moda);...
                    uni(11+(a)-moda);uni(12+(a)-moda);uni(13+(a)-moda);...
                    uni(14+(a)-moda);uni(15+(a)-moda);uni(16+(a)-moda)];
    
end

end

