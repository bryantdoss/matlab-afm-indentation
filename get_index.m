function i = get_index (array, value)
diffs = abs(array - value);
i = find(diffs == min(diffs));
i = i(1);
return