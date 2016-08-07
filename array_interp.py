# Interpolate an array when the index is not integer
# Testing of the code is needed

def array_interp(a, index):
    if sum(index<0)>0 or sum(index>len(a))>0:
        print('array_interp out of range!')
    index = index + 1E-8
    index_floor = np.array(np.floor(index), dtype=int)
    index_ceil = np.array(np.ceil(index), dtype=int)
    result = a[index_floor]*(index_ceil-index) + a[index_ceil]*(index - index_floor)
    print((index - index_floor))
    print((index_ceil - index))
    return(result)
