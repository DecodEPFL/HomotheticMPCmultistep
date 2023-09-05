function block_diag = blockDiagonal(A, n)
    block = repmat({A}, n, 1);
    block_diag = blkdiag(block{:});
end