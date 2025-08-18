
import random   
min_30_digit = 10**34
max_30_digit = (10**35) - 1 
my_list = list(range(500))

print("[", end="")
for x in range (0,499): 
    my_list[x] = random.randint(min_30_digit, max_30_digit)
    print(str(my_list[x])+",", end="")
print(str(random.randint(min_30_digit, max_30_digit))+"]")
