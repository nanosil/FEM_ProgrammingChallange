function mid = binarySearch(arr,start,en,key)
%Recursive binarySearch function. Time complexity of O(log(n)).
   if (en >= start)
        mid = floor(start + (en - start)/2);
 
        % if key is found at the middle
        if (arr(mid,1) == key)
            return;
        end
		%if key is smaller than middle element
        if (arr(mid,1) > key)
            mid = binarySearch(arr, start, mid-1, key);
            return;
        end
		%if key is larger than middle element
        mid = binarySearch(arr, mid+1, en, key);
        return;
   end
   %return here if not found
   mid =-1;
   return
   