function y=symmetricize(sol)

% Function to symmetricize solution

xp = length(sol.x);

symsol = sol;
symsol.sol(:, 1:xp, :) = sol.sol(:,:,:);

symsol.x(1:xp) = sol.x;
for ii=2:xp
symsol.sol(:, xp+ii-1, :) = symsol.sol(:, xp-ii+1, :);
symsol.x(ii+xp-1) = 2*sol.x(end) - sol.x(xp-ii+1);
end
y=symsol;

end