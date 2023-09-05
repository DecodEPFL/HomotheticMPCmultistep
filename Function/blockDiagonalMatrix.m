function result = blockDiagonalMatrix(F, p)
    % Get the size of matrix F
    [m, n] = size(F);
    
    % Initialize the result matrix with appropriate dimensions
    result = zeros(m * p, n * p);
    
    % Fill the diagonal blocks with matrix F
    for i = 1:p
        % Compute the starting indices of the current block
        startRow = (i - 1) * m + 1;
        startCol = (i - 1) * n + 1;
        
        % Copy matrix F into the current block
        result(startRow:startRow + m - 1, startCol:startCol + n - 1) = F;
    end
end
