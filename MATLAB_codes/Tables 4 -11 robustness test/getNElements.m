function [smallestNElements smallestNIdx] = getNElements(A, n)
     [ASorted AIdx] = sort(abs(A),'descend');
     smallestNElements = ASorted(1:n);
     smallestNIdx = AIdx(1:n);
end

