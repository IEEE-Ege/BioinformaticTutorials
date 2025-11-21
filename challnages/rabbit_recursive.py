def rabbit_calc(n, k):
    # Base Case: Months 1 and 2 always have 1 pair
    if  n ==1 or n ==2:
        return 1
    
    # Return (Previous Month) + (Litter Size * 2 Months Ago)
    return rabbit_calc(n-1, k) + k * rabbit_calc(n-2, k) 


print(rabbit_calc(5,3)) # k means each pair of rabbits produce 3 pairs of offspring

# if k = 1, it means each pair of rabbits produce 1 pair of offspring
