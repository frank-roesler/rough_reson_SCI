function [c4n_new, n4e_new] = delete_row(c4n,n4e,row)
    c4n_new = c4n;
    c4n_new(row,:) = [];
    delete_indices = [];
    for i=1:size(n4e,1)
        if any(n4e(i,:)==row)
            delete_indices = [delete_indices,i];
        end
    end
    keep_indices = setdiff(1:size(n4e,1), delete_indices);
    n4e_new = n4e(keep_indices,:);
    n4e_new(n4e_new>row) = n4e_new(n4e_new>row) - 1;
end