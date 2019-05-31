function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 1) pose rotation
%
% IMPORTANT NOTE: While theoretically you can use the x directly when solving
% for the P = [R t] matrix then use the K matrix to correct the error, this is
% more numeically unstable, and thus it is better to calibrate the x values
% before the computation of P then extract R and t directly

[N, ~] = size(x);
X = [X, ones(N, 1)];

% construct matrix
A = zeros(3*N, 12);
for i = 1: N
    
    Xi = X(i, :);
    ui = x(i, 1);
    vi = x(i, 2);
    
    A(i*3-2, 5:8) = -Xi;
    A(i*3-2, 9:12) = vi*Xi;
    A(i*3-1, 1:4) = Xi;
    A(i*3-1, 9:12) = -ui*Xi;
    A(i*3, 1:4) = -vi*Xi;
    A(i*3, 5:8) = ui*Xi;
    
end

[~, ~,  V] = svd(A);
P = reshape(V(:, end), [4,3])';

%RT = inv(K) * P;
RT = K \ P;
RT = RT / RT(end);

R = RT(:, 1:3);
T = RT(:, 4);

[U, S, V] = svd(R);

if det(U*V') > 0
    R = U * V';
    T = T / S(1);
else
    R = - U * V';
    T = - T / S(1);
end

C = -R' * T;

end










