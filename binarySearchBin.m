% Peter function

function [idx] = binarySearchBin(sorted_array, num)
    left = 1;
    right = length(sorted_array);
    flag = 0;

    while (right - left > 1)
        mid = ceil((left + right) / 2);

        if (num == sorted_array(mid))
            idx = mid;
            flag = 1;
            break;
        elseif (num < sorted_array(mid))
            right = mid;
        else
            left = mid;
        end
    end

    if (flag == 0)
        if (sorted_array(left) <= num && num < sorted_array(right))
            idx = left;
        else
            idx = -1;
        end
    end
end
